
# coding: utf-8

# In[1]:

import starlibrary
import glob
import pickle 


# In[ ]:

millerscalo = '345f0dd6-a9f8-11e4-9a48-90e2ba0993e2.ensemble'
sal = '8e00e316-aa91-11e4-9a48-90e2ba0993e2.ensemble'
pl215 = 'cfcb0892-ab33-11e4-9a48-90e2ba0993e2.ensemble'
pl195 = 'fe533ea8-ab8a-11e4-9a48-90e2ba0993e2.ensemble'
pl175 = '4b520882-ac4f-11e4-b168-90e2ba0993e2.ensemble'
thousking = 'a9525ffe-ac81-11e4-b168-90e2ba0993e2.ensemble' 
thousplum = 'e2c38720-ac83-11e4-b168-90e2ba0993e2.ensemble' 
quick = '12ce0e1c-ac80-11e4-b168-90e2ba0993e2.ensemble'


# In[7]:

def listclusters():
    print 'millerscalo \nsal \npl215 \npl195 \npl175 \nthousking \nthousplum \nquick' 


# In[47]:

def clusterdetails(cluster):
    ensfiles= glob.glob('./*.ensemble')
    for ensf in ensfiles:
        ensfile = open(ensf, 'rb')
        if cluster == ensf[2:]:
            runslist = pickle.load(ensfile)
            nruns = len(runslist)
            nstars = runslist[0].nstars
            if runslist[0].kingmodel:
                model = "King"
            else:
                model = "Plummer"
            print "%s -- %s model; %d runs %d stars " %(ensf, model, nruns, nstars)

        else:
            pass

    ensfile.close()

