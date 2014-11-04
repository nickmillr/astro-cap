# File: starlibrary.py (modified from starlab.py)
# Name: Nicholas Miller
# Date: 10.10.14
# Desc: Program provides the necessary tools for cluster simulations
# 		(modified form Dr. Cavendish McKay's starlab.py file
# Usage: The program reads input containing details for cluster structure.
#        The program prompts cpus to simulate clusters and export output.

# <Imports>

from subprocess import Popen, PIPE
import os
import datetime
import re
import numpy as np
import pandas as pd
from lxml import etree
import uuid
import glob
import pickle
#import pylab as plt


warehousepath = ".warehouse/"   


# <Class creation>

class Story:
    def __init__(self):
        self.story_lines = []
        self.story_vals = dict()
        self.story_subobjects = []
        self.kind = None
        return
    def __str__(self):
        return "%s, %d lines, %d values, %d subobjects" % (self.kind, len(self.story_lines), len(self.story_vals.keys()), len(self.story_subobjects))
    def process_line(self, line):
        # if we can tokenize into an equality, store in the dict, otherwise as a line.
        chunks = re.split('=', line)
#        print chunks, len(chunks)
        if len(chunks) == 2:
            self.story_vals[chunks[0].strip()] = chunks[1].strip()
        else:
            self.story_lines.append(line)

# <set up the Run>
class Run:
    def __init__(self, kingmodel=True, w0=1.5, nstars=2500, masstype=1, runlength=100, exponent=-2.35, binarypercent=.1, binarypoplower = 1.0, binarypopupper = 10):

        self.kingmodel = kingmodel
        self.w0 = w0
        self.nstars = nstars

        # mass scaling
        self.masstype=masstype
        self.exponent=exponent
        self.lowerlimit = 0.1
        self.upperlimit = 20

        # binary scaling
        self.binarypercent = binarypercent
        self.binarylimit = 0.25
        self.binarypoplower = binarypoplower
        self.binarypopupper = binarypopupper

        # kira parameters
        self.runlength=runlength
        self.diagout = 0.5
        self.dumpout = 0.5
        
# <codecell>
def parse_output(results):
    stories = []
    lines = re.split("\n", results)

    nextidx = -1
    
    for index, line in enumerate(lines):
        if index >= nextidx:
            storystart = re.match("^\((\w+)",line)
            if storystart:
                #print "in parse_output, calling parse_lines:", index, storystart.group(1)
                nextidx, newstory = parse_lines(lines, index+1, storystart.group(1))
                #print "in parse_output, back from parse_lines:", nextidx, str(newstory)
                stories.append(newstory)
    return stories

# <codecell>

def parse_lines(lines, startidx, name):
    thestory = Story()
    thestory.kind = name
    nextidx = -1
    for index,line in enumerate(lines[startidx:]):
        if index >= nextidx-startidx:
            #print "%s Line is %s"%(name, line)
            storystart = re.match("^\((\w+)",line)
            storyend = re.match("\)%s"%name, line)
            if storyend: # we've hit the end of our story; get out and pass it back up
                endindex = index
                break
            elif storystart: # new story; start up a new parse_lines
                #print "in parse_lines, calling parse_lines:", startidx+index+1, storystart.group(1)
                nextidx, newstory = parse_lines(lines, startidx+index + 1, storystart.group(1))
                #print "back", nextidx, str(newstory)
                thestory.story_subobjects.append(newstory)
            else:
                thestory.process_line(line)
    return endindex+startidx+1, thestory

# <codecell>

def extract_particle_dynamics(story):
    """ recursively extract dynamics objects to get masses, positions, and velocities.

    Be careful of particles with subparticles."""
    particles = []
    if story.kind == "Particle":
        if len(story.story_subobjects) > 4:
            # get r, v, as center of mass quantities
            for subobj in story.story_subobjects:
                if subobj.kind == 'Dynamics':
                    xcom, ycom, zcom = subobj.story_vals['r'].split(" ")
                    vxcom, vycom, vzcom = subobj.story_vals['v'].split(" ")
            # get relative positions and velocities of sub-particles
            subparticles = []
            for subobj in story.story_subobjects:
                if subobj.kind == 'Particle':
                    subparticles.extend(extract_particle_dynamics(subobj))
                    
            # add COM r, v to sub-particles, and append
            for particle in subparticles:
                particles.append((particle[0] + float(xcom),
                                  particle[1] + float(ycom),
                                  particle[2] + float(zcom),
                                  particle[3] + float(vxcom),
                                  particle[4] + float(vycom),
                                  particle[5] + float(vzcom),
                                  particle[6]))
        else: # only 4 subobjects, so this is an individual star
            for subobj in story.story_subobjects:
                if subobj.kind == 'Dynamics':
                    x,y,z = subobj.story_vals['r'].split(" ")
                    vx,vy,vz = subobj.story_vals['v'].split(" ")
                    m = subobj.story_vals['m']
                    particles.append((float(x), float(y), float(z), float(vx), float(vy), float(vz), float(m)) )
    return particles
    
# <codecell>

def vis_story_3d(story_list):
    """visualize a story list. """
    
    xvals = []
    yvals = []
    zvals = []
    vxs = []
    vys = []
    vzs = []
    masses = []
    
    for story in story_list:
        partlist = extract_particle_dynamics(story)
        
    for particle in partlist:
        xvals.append(particle[0])
        yvals.append(particle[1])
        zvals.append(particle[2])
        vxs.append(particle[3])
        vys.append(particle[4])
        vzs.append(particle[5])
        masses.append(particle[6])
        
    # now do the plot
    fig = plt.figure(figsize=(10,10))
    ax = fig.gca(projection='3d')
    ax.plot(np.array(xvals), np.array(yvals), np.array(zvals), ".")

# <codecell>

def premain(startn):
    """Run a plummer model for 10 dynamical times and return the number of stars remaining."""
    from subprocess import Popen, PIPE
    from starlibrary import parse_output, extract_particle_dynamics
    import random
    
    seed = random.randint(0,9999999999)
    
    print "running %d particles" % startn
    cmds = []

    cmds.append(["makeking", "-n", "%d"%startn, "-w", "5", "-i",  "-u", "-s", "%d"%seed])
    cmds.append(["makemass", "-f", "2", "-l", "0.1,", "-u", "20"])
    cmds.append(["makesecondary", "-f", "0.1", "-l", "0.25"])
    cmds.append(["makebinary", "-l", "1", "-u", "10"])
    cmds.append(["scale", "-m", "1", "-e", "-0.25", "-q", "0.5"]) 
    cmds.append(["kira", "-t", "100", "-d", "1", "-D", "2", "-f", "0.3", "-n", "10", "-q", "0.5", "-G", "2", "-B"])

    procs = []
    for index, cmd in enumerate(cmds):
        print index, cmd
        if index > 0:
            procs.append(Popen(cmd, stdout=PIPE, stderr=PIPE, stdin=procs[index-1].stdout))
        else:
            procs.append(Popen(cmd, stdout=PIPE, stderr=PIPE))
    inp = procs[-1].stdout
    
    result = procs[-1].communicate()
    slist = parse_output(result[0])
    return len(extract_particle_dynamics(slist[-1]))

# <codecell>

def grab_energies(slist):
    times = []
    energies = []
    kinetics = []
    potentials = []
    for story in slist:
        times.append(float(slist.story_subobjects[1].story_vals['t']))
        energies.append(float(slist.story_subobjects[1].story_vals['total_energy']))
        kinetics.append(float(slist.story_subobjects[1].story_vals['kinetic_energy']))
        potentials.append(float(slist.story_subobjects[1].story_vals['potential_energy']))
    return times, energies, kinetics, potentials
# <codecell>

def l2norm(E, U, T):
    U0 = 2*E[0]
    T0 = -E[0]
    l2U = sum((np.array(U)-U0)**2)
    l2T = sum((np.array(T)-T0)**2)
    return l2U, l2T

# <codecell>

def fs_process_commands(cmds):
    """Process commands, storing the results on the filesystem."""
    procs = []
    for index, cmd in enumerate(cmds):
        if index > 0:
            procs.append(Popen(cmd, stdout=PIPE, stdin=procs[index-1].stdout))
        else:
            procs.append(Popen(cmd, stdout=PIPE))
        inp = procs[-1].stdout
    result = procs[-1].communicate()

    theuuid = uuid.uuid1()
    storename = str(theuuid)+".kiraout"

    #if not os.path.isdir(warehousepath):
    #    os.makedirs(warehousepath)

    store = open(warehousepath+storename, "w")
        
    store.write(result[0])
    store.close()
    
    return theuuid

# <codecell>
def fs_panel(storename):
    """Grab a story list from the filesystem and turn it into a pandas panel."""
    infile = open(storename, 'rb')
    partstart = re.compile("^\(Particle")
    partend = re.compile("^\)Particle")
    dynstart = re.compile("^\(Dynamics")
    dynend = re.compile("^\)Dynamics")

    snap = 0
    insnap = False
    panbase = pd.Panel()
    
    for line in infile:
        if partstart.match(line):
            pass
        elif partend.match(line):
            pass
        elif dynstart.match(line):
            pass
        elif dynend.match(line):
            pass
        
    return pan


# <codecell>
def kiraout_to_xml(theuuid):
    storebase = str(theuuid)
    storename = str(theuuid)+".kiraout"
    
    xmlframes = []

    frame = 0
    nframes = 0
    level = 0
    
    store = open(warehousepath+storename, "rb")

    storystart = re.compile("^\((\w+)")
    storyend = re.compile("^\)(\w+)")
    illegalstart = re.compile("^-")

    for line in store:
        if storystart.match(line):
            #print "entering: %s, %d" % (storystart.match(line).group(1), level)
            if level==0:
                nframes += 1
                xmlname = storebase + ".%d.xml" % frame
                xmlframes.append(xmlname)
                xmlfile = open(warehousepath+xmlname, "wb")
                xmlfile.write('<?xml version="1.0" encoding="ISO-8859-1"?>')
            level += 1
            xmlfile.write("<%s>\n"%storystart.match(line).group(1))
        elif storyend.match(line):
            #print "leaving: %s, %d" % (storyend.match(line).group(1), level)
            level -= 1
            xmlfile.write("</%s>\n"%storyend.match(line).group(1))
            if level == 0:
                frame += 1
                xmlfile.close()
        else:
            chunks = re.split("=", line)
            if len(chunks) == 2:
		# xmlfile.write("<value %s=\"%s\" />\n" %(re.sub("\,* ","_",chunks[0].strip()), chunks[1].strip()))
                # xml attributes can't begin with punctuation or a digit. Make sure we abide by that rule.
		if re.match("^\s*[a-zA-Z]", chunks[0]):
		     xmlfile.write("<value %s=\"%s\" />\n" %(re.sub("\,* ","_",chunks[0].strip()), chunks[1].strip()))
		else:
                     xmlfile.write(line)
            else:
                xmlfile.write(line)

    store.close()
    return nframes

# <codecell>


def quantile_radius(pan, quant):
    """Get the average radius for a given mass quantile.

    Compares the average radius for the given mass quantile 0% -> quant
    with that of the complement ((1-quant) -> 100%)

    Takes inputs:
    pan -- pandas panel, where the item index is sorted by mass
    quant -- the quantile, in percent, to use.

    returns a tuple of time series:
    (lightradius, heavyradius, light/heavy) 
    """
    indices = sorted(pan.keys(), key=int)
    number = int(quant*len(indices))
    print "%d runs in requested quantile" % number
    heavyslice = [pan[index] for index in indices[:number]]
    lightslice = [pan[index] for index in indices[-number:]]
    
    print "heavy masses:", [star['m'][0] for star in heavyslice]
    print "light masses:", [star['m'][0] for star in lightslice]
    
    heavyr = np.zeros_like(heavyslice[0]['radius'].values)
    for star in heavyslice:
        heavyr += star['radius'].values
    heavyr /= number
    
    lightr = np.zeros_like(lightslice[0]['radius'].values)
    for star in lightslice:
        lightr += star['radius'].values
    lightr /= number
    
    ratio = lightr/heavyr
    return lightr, heavyr, ratio

# <codecell>

from matplotlib import animation
def animate_panel_old(pan, prng=10.0, filebase='basic_animation'):
    """Turn a panel of star cluster evolution snapshots into a 3d animation.

    Arguments are:
    pan: the panel
    prng: plot range for all three axes.
    filebase: the base name of the file for the animation.

    returns:
    Nothing."""
    
    fig = plt.figure(figsize=(10,10))
    fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
    ax = fig.gca(projection='3d')
    ax.set_autoscale_on(False)
    
    nframes = pan.shape[1]
    starfield, = ax.plot([], [], [],'bo', ms=7)
    ax.set_xlim3d(-prng, prng)
    ax.set_ylim3d(-prng, prng)
    ax.set_zlim3d(-prng, prng)

    def init():
        """initialize animation"""
        starfield.set_data([], [])
        starfield.set_3d_properties(zs=[])
        return starfield
    
    def doframe(n):
        xvals = []
        yvals = []
        zvals = []
        masses = []
        
        for starid in pan.keys():
            xvals.append(pan[starid]['r'][n][0])
            yvals.append(pan[starid]['r'][n][1])
            zvals.append(pan[starid]['r'][n][2])
            masses.append(pan[starid]['m'][n])
    # not super happy with the colors yet.
        starfield.set_data(np.array(xvals), np.array(yvals))
        starfield.set_3d_properties(zs=np.array(zvals))
        return starfield

    anim = animation.FuncAnimation(fig, doframe,
                               frames=nframes, interval=20, blit=True)
    anim.save("%s.mp4"%filebase, fps=30)
    #ax.plot(np.array(xvals), np.array(yvals), np.array(zvals), ".")


def animate_from_fs(filebase, nframes, prng=10.0, use_warehouse=True):
    import pylab as plt
    from mpl_toolkits.mplot3d import axes3d
    fig = plt.figure(figsize=(10,10))
    fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
    ax = fig.gca(projection='3d')
    ax.set_autoscale_on(False)

    mcomfield, = ax.plot([], [], [],'go', ms=7)
    ax.set_xlim3d(-prng, prng)
    ax.set_ylim3d(-prng, prng)
    ax.set_zlim3d(-prng, prng)

    def init():
        """initialize animation"""
        mcomfield.set_data([], [])
        mcomfield.set_3d_properties(zs=[])
        return (mcomfield)
    
    def doframe(n):
        xvals = []
        yvals = []
        zvals = []
        masses = []
        
        if use_warehouse:
            xmlname = warehousepath + filebase + ".%d.xml"% n
        else:
            xmlname = filebase + ".%d.xml"% n
            
        dataframe, time, nparticles = process_frame(xmlname)
        
        xvals = dataframe.xs('x').values
        yvals = dataframe.xs('y').values
        zvals = dataframe.xs('z').values
        
        mcomfield.set_data(dataframe.xs('x_mcom').values, dataframe.xs('y_mcom').values)
        mcomfield.set_3d_properties(zs=np.array(dataframe.xs('z_mcom').values))
        return (mcomfield)

    anim = animation.FuncAnimation(fig, doframe,
                               frames=nframes, interval=20, blit=True)
    moviename = "%s.mp4"%filebase
    anim.save(moviename, fps=30, bitrate=4000)
    return moviename

def str_to_vect(vectstr):
    return np.array(vectstr.split(" "), float)

def process_frame(xmlname):
    frame = etree.parse(xmlname)
    nparticles = int(frame.xpath('/Particle/value[@ N]/@ N')[0])
    time = float(frame.xpath('/Particle/Dynamics/value[@ t]/@ t')[0])
    
    center_r = str_to_vect(frame.xpath('/Particle/Dynamics/value[@ r]/@ r')[0])
    com_pos = str_to_vect(frame.xpath('/Particle/Dynamics/value[@ com_pos]/@ com_pos')[0])
    mcom_pos = str_to_vect(frame.xpath('/Particle/Dynamics/value[@ mcom_pos]/@ mcom_pos')[0])
    mass_scale = float(frame.xpath('/Particle/Star/value[@ mass_scale]/@ mass_scale')[0])
    time_scale = float(frame.xpath('/Particle/Star/value[@ time_scale]/@ time_scale')[0])
    size_scale = float(frame.xpath('/Particle/Star/value[@ size_scale]/@ size_scale')[0])

    particledict = {}
    for part in frame.xpath('/Particle/Particle'):
        subpartdict = flatten_particle(part)
        particledict.update(subpartdict)
    for key in particledict.keys():
        particledict[key]['x_mcom'] = particledict[key]['x'] - mcom_pos[0] + center_r[0]
        particledict[key]['y_mcom'] = particledict[key]['y'] - mcom_pos[1] + center_r[1]
        particledict[key]['z_mcom'] = particledict[key]['z'] - mcom_pos[2] + center_r[2]
        particledict[key]['x_com'] = particledict[key]['x'] - com_pos[0] + center_r[0]
        particledict[key]['y_com'] = particledict[key]['y'] - com_pos[1] + center_r[1]
        particledict[key]['z_com'] = particledict[key]['z'] - com_pos[2] + center_r[2]
    framedataframe = pd.DataFrame(particledict)
    
    return framedataframe, time, nparticles

def xml_to_panel(basename, nframes):
    paneldict = {}
    for frame in xrange(nframes):
        xmlname = basename + ".%d.xml"%frame
        df, time, nparticles = process_frame(xmlname)
        paneldict[time] = df
    thepanel = pd.Panel(paneldict)
    return thepanel.transpose(items='minor', major='items', minor='major')
    



def process_commands(cmds):
    procs = []
    for index, cmd in enumerate(cmds):
        if index > 0:
            procs.append(Popen(cmd, stdout=PIPE, stderr=PIPE, stdin=procs[index-1].stdout))
        else:
            procs.append(Popen(cmd, stdout=PIPE, stderr=PIPE))
        inp = procs[-1].stdout
    
    result = procs[-1].communicate()
    return result

def story_list_to_panel(slist):
    nsnaps = len(slist)
    nstars = slist[0].story_vals['N']
    dict_for_panel = {}
    
    for snapnum, snapstory in enumerate(slist):
        partlist = flatten_parts(snapstory.story_subobjects[4:], snapstory.story_subobjects[1])
        for part in partlist:
            # make the star dict if necessary
            starlabel = part["i"]
            pos = part['r']
            radius = np.sqrt(pos.dot(pos))
            if dict_for_panel.get(starlabel,None) == None:
                dict_for_panel[starlabel] = {}
                for key in part.keys():
                    dict_for_panel[starlabel][key] = [ part[key] ]
                dict_for_panel[starlabel]['radius'] = [ radius ]
            else:
                for key in part.keys():
                    try:
                        dict_for_panel[starlabel][key].append(part[key])
                    except KeyError:
                        dict_for_panel[starlabel][key] = [ None ] * snapnum
                        dict_for_panel[starlabel][key].append(part[key])
                dict_for_panel[starlabel]['radius'].append(radius)
            # check to see if anything is missing...
            reflength = len(dict_for_panel[starlabel]['radius'])
            for key in dict_for_panel[starlabel].keys():
                if len(dict_for_panel[starlabel][key]) < reflength:
                       #print "Whoa. problem. i=%s, t=%s, key=%s" % (starlabel, snapnum, key)
                       dict_for_panel[starlabel][key].append(None)

#    print "number of stars: ",len(dict_for_panel.keys())
#    for key in dict_for_panel.keys():
#        lengthlist = [len(dict_for_panel[key][subkey]) for subkey in dict_for_panel[key].keys()]
#        print key, lengthlist
#        print dict_for_panel[key].keys()
#        print "-------------------------\n\n"
        
    thepanel = pd.Panel(dict_for_panel)
    return thepanel

def flatten_parts(partlist, parentinfo):
    flattened = []
    for partnum, part in enumerate(partlist):
        if len(part.story_subobjects) > 4:
            flattened.extend(flatten_parts(part.story_subobjects[4:], part.story_subobjects[1]))
        else:
            partdict = dict(part.story_vals, **part.story_subobjects[1].story_vals)
            flattened.append(partdict)
    # deal with positions
    pos = np.array(parentinfo.story_vals['r'].split(" "), float)
    for part in flattened:
        try:
            part['r'] = part['r'] + pos
        except TypeError:
            part['r'] = np.array(part['r'].split(" "), float) + pos
    return flattened


def single_run(therun):
    """Perform a single run from a Run object."""
    cmds = []

    if therun.kingmodel:
        cmds.append(["makeking",
                      "-n", "%d"%therun.nstars,
                      "-w", "%3.1f"%therun.w0,
                      "-i"])
    else:
        cmds.append(["makeplummer",
                      "-n", "%d"%therun.nstars,
                      "-i"])

    cmds.append(["makemass",
                  "-f", "%d"%therun.masstype,
                  "-l", "%f"%therun.lowerlimit,
                  "-u", "%f"%therun.upperlimit,
                  "-i"])

    cmds.append(["makesecondary",
                  "-f", "%d"%therun.binarypercent,
                  "-l", "%f"%therun.binarylimit
		])

    cmds.append(["scale",
                  "-m", "1",
                  "-e", "-0.25",
                  "-q", "0.5"])

    cmds.append(["makebinary",
                  "-l", "%d"%therun.binarypoplower,
                  "-u", "%f"%therun.binarypopupper
		])

    cmds.append(["kira", "-t", "%d"%therun.runlength,
                  "-d", "%f"%therun.diagout,
                  "-D", "%f"%therun.dumpout,
                  "-n", "10",
                  "-q", "0.5"])

    therun.uuid = fs_process_commands(cmds)
    therun.nframes = kiraout_to_xml(therun.uuid)
    return therun

def animate_run(therun):
    theuuid = therun.uuid
    moviename = animate_from_fs(str(theuuid), therun.nframes, use_warehouse=True)
    return moviename

def flatten_particle(part):
    """Turn a prarticle into a dictionary with the appropriate values, recursing if necessary."""
    nparticles = int(part.xpath('./value[@ N]/@ N')[0])
    #print nparticles
    r = str_to_vect(part.xpath('./Dynamics/value[@ r]/@ r')[0])
    v = str_to_vect(part.xpath('./Dynamics/value[@ v]/@ v')[0])
    
    if nparticles == 1:
        parti = int(part.xpath('./value[@ i]/@ i')[0])
        m = float(part.xpath('./Dynamics/value[@ m]/@ m')[0])
        radius = np.sqrt(r.dot(r))
        pot = float(part.xpath('./Dynamics/value[@ pot]/@ pot')[0])
        partdict = {parti: {'x':r[0], 'y':r[1], 'z':r[2], 'm':m, 'radius':radius} }
        partdict[parti].update({'vx':v[0], 'vy':v[1], 'vz':v[2], 'pot':pot})
    else:
        # this node is not a leaf; we have to do the subtree
        partdict = {}
        for subpart in part.xpath('./Particle'):
            tpd = flatten_particle(subpart)
            #print tpd
            partdict.update(tpd)
        for key in partdict.keys():
            partdict[key]['x'] += r[0]
            partdict[key]['y'] += r[1]
            partdict[key]['z'] += r[2]
            partdict[key]['vx'] += v[0]
            partdict[key]['vy'] += v[1]
            partdict[key]['vz'] += v[2]
            partdict[key]['radius'] = np.sqrt(partdict[key]['x']**2+partdict[key]['y']**2+partdict[key]['z']**2)
            partdict[key]['kinetic'] = (partdict[key]['vx']**2+partdict[key]['vy']**2+partdict[key]['vz']**2)*partdict[key]['m']/2.0
    return partdict

def extract_from_frame(therun, framenum):
    filebase = str(therun.uuid)
    xmlname = '.warehouse/' + filebase + ".%d.xml" % framenum
    frame = etree.parse(xmlname)

    extracted = {}
    
    nparticles = int(frame.xpath('/Particle/value[@ N]/@ N')[0])
    time = float(frame.xpath('/Particle/Dynamics/value[@ t]/@ t')[0])

    center_r = str_to_vect(frame.xpath('/Particle/Dynamics/value[@ r]/@ r')[0])
    com_pos = str_to_vect(frame.xpath('/Particle/Dynamics/value[@ com_pos]/@ com_pos')[0])
    mcom_pos = str_to_vect(frame.xpath('/Particle/Dynamics/value[@ mcom_pos]/@ mcom_pos')[0])
    mass_scale = float(frame.xpath('/Particle/Star/value[@ mass_scale]/@ mass_scale')[0])
    time_scale = float(frame.xpath('/Particle/Star/value[@ time_scale]/@ time_scale')[0])
    size_scale = float(frame.xpath('/Particle/Star/value[@ size_scale]/@ size_scale')[0])

    particledict = {}
    for part in frame.xpath('/Particle/Particle'):
        subpartdict = flatten_particle(part)
        particledict.update(subpartdict)
    for key in particledict.keys():
        particledict[key]['x_mcom'] = particledict[key]['x'] - mcom_pos[0] + center_r[0]
        particledict[key]['y_mcom'] = particledict[key]['y'] - mcom_pos[1] + center_r[1]
        particledict[key]['z_mcom'] = particledict[key]['z'] - mcom_pos[2] + center_r[2]
        particledict[key]['x_com'] = particledict[key]['x'] - com_pos[0] + center_r[0]
        particledict[key]['y_com'] = particledict[key]['y'] - com_pos[1] + center_r[1]
        particledict[key]['z_com'] = particledict[key]['z'] - com_pos[2] + center_r[2]
        particledict[key]['radius'] = np.sqrt(particledict[key]['x_mcom']**2 +
                                              particledict[key]['y_mcom']**2 +
                                              particledict[key]['z_mcom']**2)
    framedataframe = pd.DataFrame(particledict)

    return framedataframe, time, nparticles

def run_to_panel(therun):
    nframes = int(therun.runlength/therun.dumpout)
    paneldict={}
    for frame in xrange(nframes):
        df, time, nparticles = extract_from_frame(therun, frame)
        paneldict[time] = df
    thepanel = pd.Panel(paneldict)
    return thepanel.transpose(items='major', major='minor', minor='items')

def massradius(pan, quantity='radius', fraction=0.1):
    times = pan[quantity].keys()
    high_mass_radius_list = []
    low_mass_radius_list = []
    nstars = pan[quantity][times[0]].count()
    nfrac = int(nstars * fraction)
    
    hmseries = pan[quantity][0:nfrac].median()
    lmseries = pan[quantity][-nfrac:-1].median()
    ratioseries = lmseries/hmseries
        #print time, hmr, lmr
        
#    lmseries = pd.Series(data=low_mass_radius_list, index=times)
#    hmseries = pd.Series(data=high_mass_radius_list, index=times)
#    ratioseries = pd.Series(data=ratios, index=times)
    return(ratioseries, lmseries, hmseries)

def list_ensembles():
    ensfiles= glob.glob('./*.ensemble')
    print "%d ensembles:\n" % len(ensfiles)
    for ensf in ensfiles:
        ensfile = open(ensf, 'rb')
        try:
            runslist = pickle.load(ensfile)
            nruns = len(runslist)
            nstars = runslist[0].nstars
            if runslist[0].kingmodel:
                model = "King"
            else:
                model = "Plummer"
            print "%s -- %s model; %d runs %d stars " %(ensf, model, nruns, nstars)
        except EOFError:
            print "%s -- No data" %ensf
        ensfile.close()


def avg_series(ensname, output=False):
    storename = "%s-avg.h5"%ensname
    store = pd.HDFStore("%s-avg.h5"%ensname)
    try:
        ratiodf = store['ratiodf']
        lowdf = store['lowdf']
        highdf = store['highdf']
    except KeyError:
        ensfile = open(ensname, 'rb')
        runlist = pickle.load(ensfile)

        rsdict = {}
        lsdict = {}
        hsdict = {}
        i = 0
        for run in runlist:
            key = "%d"%i
            if output:
                print key
            try:
                thepan = run_to_panel(run)
                rseries, lmseries, hmseries = massradius(thepan)
                rsdict[key] = rseries
                lsdict[key] = lmseries
                hsdict[key] = hmseries
            except IOError:
                pass
            i+=1
        highdf = pd.DataFrame(hsdict)
        lowdf = pd.DataFrame(lsdict)
        ratiodf = pd.DataFrame(rsdict)
    
        # compute averages
        highdf['avg'] = highdf.mean(axis=1)
        lowdf['avg'] = lowdf.mean(axis=1)
        ratiodf['avg'] = ratiodf.mean(axis=1)
    
        # compute stddev
        highdf['stddev'] = highdf.std(axis=1)
        lowdf['stddev'] = lowdf.std(axis=1)
        ratiodf['stddev'] = ratiodf.std(axis=1)
 
        store['highdf'] = highdf
        store['lowdf'] = lowdf
        store['ratiodf'] = ratiodf
        
    return ratiodf, lowdf, highdf

def r_and_m(clusterstep):
    rvals = []
    masses = []
    
# determine the radii and masses of stars for a timestep in kira output
# uses extract_particle_dynamics to open story and create lists of stars
    
    partlist = extract_particle_dynamics(clusterstep)
    
    for particle in partlist:
        masses.append(particle[6])
        x = particle[0]
        y = particle[1]
        z = particle[2]
        r = np.sqrt(x*x + y*y + z*z)
        rvals.append(r)
#        print r, particle[6]
        
    return rvals, masses

def medianpos(mlow, mhigh, rad, masses):
    index = 0
    r_lowmass = []
    r_highmass = []
    
# find median radius for stars below a low mass cut-off and above
# a high-mass cut-off given the desired cut-offs as well as lists of
# radii and masses (ordered the same)
    
    for mass in masses:
        if mass <= mlow:
            r_lowmass.append(rad[index])
        if mass >= mhigh:
            r_highmass.append(rad[index])
        index = index + 1
        
    medlow = np.median(r_lowmass)
    medhigh = np.median(r_highmass)
    
    return medlow, medhigh

def medianpos_percentile(storylist):
    masscut = 0.1
    first = 1
    median_low = []
    median_high = []
    
# takes a storylist generated from kira output and produces a list of
# median star positions for stars in the upper & lower percentile as
# set in the parameter "masscut", above

# when looking at the initial timestep, it uses the percentile to figure
# out actual cut-offs in mass.  This way, we do not have to worry about 
# evaporation changing which stars fall in our observed bins.  
    
    for story in storylist:
        (radii, masses) = r_and_m(story)
        if first == 1:
            masssort = np.sort(masses)
            nummass = len(masses)
            indexlo = int(np.floor(masscut * nummass))
            indexhi = int(np.floor(nummass - indexlo))
            masslow = masssort[indexlo]
            masshigh = masssort[indexhi]
            print masslow, masshigh
            first = 0
        (medlow, medhigh) = medianpos(masslow, masshigh, radii, masses)
        median_low.append(medlow)
        median_high.append(medhigh)
    return median_low, median_high
    
