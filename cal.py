import pandas, numpy,scipy,os
import datetime, dateparser
import pylab
pylab.ion() 
import pygecros_fen as pyg
#from pycrop import gecros_utils
import seaborn

    
def read_wf_csiro(folder='../AgMIP_Cal_phase3/CSIRO_data_set/'):
    data={}
    files=[[x,x[:-8]] for x in os.listdir(folder) if x[-4:]=='.met']
    for i in files:
        fname=folder+i[0]
        print (fname)
        f=open(fname,'r')
        l=f.readlines()
        istart=[i for i in range(len(l)) if '!' in l[i] or 'year' in l[i]]
        #latitude=[i for i in range(len(l)) if 'latitude' in l[i] ]
        f.close()
        
        if i[1] not in data.keys():
            data[i[1]]=pandas.DataFrame()
        names=pandas.read_table(fname, delimiter=' ', skipinitialspace=True,skiprows=istart[-1],nrows=0).columns 
        w=pandas.read_table(fname, delimiter=' ',skipinitialspace=True, skiprows=istart[-1]+2,names=names)  
        w.index=pandas.to_datetime(w.apply(lambda x: datetime.date(int(x.year),1,1)+datetime.timedelta(days=int(x.day)-1),axis=1))
        data[i[1]]=pandas.concat([data[i[1]],w])
        print ('ok')
    for i in data.keys():data[i].sort_index(inplace=True)  
    return data
            #year  day  radn  maxt  mint  rain   pan    vp    code
#2011-01-01  2011    1  28.0  41.5  23.5   0.0  11.6  13.0  222022
#2011-01-02  2011    2  18.0  40.5  25.5   0.0   5.8  18.0  222022
class sim():
  def __init__(self,):
    self.latitude={'Eradu':-28.69,'LakeBolac':-37.71,'Minnipa':-32.84,'SpringRidge':-31.39,'Turretfield':-34.55,'Walpeup':-35.12,'Temora':-34.41,'Nangwee':-27.66,'Corrigin':-32.33,'Bungunya':-28.43}
    self.wf=read_wf_csiro()
    self.d=pandas.io.excel.read_excel('../AgMIP_Cal_phase3/CSIRO_data_set/training Zadoks dates.xlsx',sheet_name='Hoja1',header=0,index_col=[0,1,2])      

  def obs(self,):
    cal_sites=self.d[numpy.logical_and(numpy.logical_not(pandas.isna(self.d.Zadok65)),numpy.logical_not(pandas.isna(self.d.Zadok89)))]
    data={}
    for i in cal_sites.iterrows():
      year=i[0][0]
      site=i[0][1]
      toss=i[0][2]
      pld=i[1].Sowing_date
      obsZadok65=i[1].Zadok65.toordinal()
      obsZadok89=i[1].Zadok89.toordinal() 
      sit_name=str(site)+'-'+str(year)+'-'+str(toss)
      data[sit_name]=pandas.DataFrame(data={'Date':[pld.strftime("%Y-%m-%d")] , 'Zadok65':[obsZadok65] ,'Zadok89':[obsZadok89]},index=[1])
    return data
  def run_sites(self,sit_names,*args):
    data={}
    #print(sit_names,args)
    MTDV,MTDR, PSEN,VSEN,VDSA,*rest= list(args[0]) + [None]*10#.tolist()
    #print(sit_names,MTDV , MTDR, PSEN, VSEN, VDSA,rest)
    if MTDV==None:MTDV=32.
    if MTDR==None:MTDR=31.
    if PSEN==None:PSEN=0.06
    if VSEN==None:VSEN=0.003
    if VDSA==None:VDSA=50
    
    for j in sit_names:
      site,year,toss=j.split('-')
      i=self.d.query('Year=={0} and Site=="{1}" and Time_of_Sowing=="{2}"'.format(int(year), site, toss))   
      #year=i[0][0]
      #site=i[0][1]
      #toss=i[0][2]
      SLP = -1.    # SLP    = 1. for short-day crops;    = -1. for long-day crops.
      TBD=0.       # Base temperature for phenology (Tb)                                                  C
      TOD=25.0     # Optimum temperature for phenology (To)                                               C
      TCD=37.      # Ceiling temperature for phenology (Tc)                                               C
      TSEN=1.5     #Curvature for temperature response (ct)   (value range 0.5-3.0)                      --
      SPSP=0.2     #Development stage (DS) for start of fotoperiod sensitive phase (v1) (0-0.5)          --
      EPSP=0.7     #Development stage (DS) for end of photoperiod sensitive phase  (v2) (0.5-0.8)        --
      SVSP=0.      #Development stage (DS) for start of vernalization sensitive phase (v1) (0-0.5)          --
      EVSP=0.5     #Development stage (DS) for end of vernalization sensitive phase  (v2) (0.5-0.8)        --    
      INSP=-2.     # Inclination of sun angle for computing daylenght for photoperiodism(alpha)           degrees
      #MTDV=32      #Minimum thermal days for vegetative growth phase                                     d
      #MTDR=31.     #18,14#Minimum thermal days for reproductive (seed fill)phase                               d
      #PSEN=-0.06   #Photoperiod sensitivity of phenological development                                  h-1
      #VSEN= 0.003  #Vernalization sensitivity   (0.033 WW ,0.003 SWhigh latitude , 0 SW)                     d
      #VDSA=50      #winterW USA #Vernalization saturation day units                                       d
      PSEN*=-1
      pld=i.iloc[0].Sowing_date 
      STTIME= pld + datetime.timedelta(days=11)                                                                        
      FINTIM =STTIME+datetime.timedelta(days=360)    #time to finish run       
      LAT=self.latitude[site]
      WEATHER=self.wf[site]
      ET,AT,MT=pyg.run(SLP=SLP,TBD=TBD,TOD=TOD,TCD=TCD,TSEN=TSEN,SPSP=SPSP,EPSP=EPSP,SVSP=SVSP,EVSP=EVSP,INSP=INSP,MTDV=MTDV,MTDR=MTDR,PSEN=PSEN,VSEN=VSEN,VDSA=VDSA,COEFT=0.,STTIME=STTIME,FINTIM=FINTIM,WEATHER=WEATHER,LAT=LAT)
      #print (year,site,toss,ET,AT,MT)
      
      data[j]=pandas.DataFrame(data={'Date':[pld.strftime("%Y-%m-%d")] , 'Zadok65':[AT] ,'Zadok89':[MT]},index=[1])#strftime("%Y-%m-%d")
    return data
  def plot_obs(self,d):
    sit_names=d.keys()
    df=pandas.DataFrame()
    for j in sit_names:
      site,year,toss=j.split('-')
      #self.wf[site]
      df[j]=self.wf[site].query('index>"{0}" and index<"{1}"'.format(d[j].Date[1],datetime.datetime.fromordinal(d[j].Zadok89[1]).date().isoformat()) ).mean()  
    df=df.T.set_index(numpy.array([i.split('-') for i in df.columns]).T.tolist())   
    df.index.set_names(['site','year','toss'],inplace=True)    
    fig,ax=pylab.subplots(1,2)  
    seaborn.boxplot(y='mint',x=df.index.get_level_values(0),data=df,ax=ax[0])  
    seaborn.boxplot(y='maxt',x=df.index.get_level_values(0),data=df,ax=ax[1])
    return df
    
  def run_optim(self,params=[]):
    #if type(params)!=list:params=list(params)
    sit_names=[]
    #obs_emerg=[]
    #obsZadok65=[]
    #obsZadok89=[]
    cal_sites=self.d[numpy.logical_and(numpy.logical_not(pandas.isna(self.d.Zadok65)),numpy.logical_not(pandas.isna(self.d.Zadok89)))]#.loc[:,['Zadok65','Zadok89']]    #los que estan completos en ambos
    #cal_sites_all=self.d.query('(Year==2010 or Year==2011) and (Site=="Eradu" or Site=="LakeBolac" or Site=="Minnipa" or Site=="SpringRidge")').loc[:,['Zadok65','Zadok89']] #todos
    for i in cal_sites.iterrows():
      year=i[0][0]
      site=i[0][1]
      toss=i[0][2]
      sit_names.append(str(site)+'-'+str(year)+'-'+str(toss))
      #obsZadok65.append(i[1].Zadok65.toordinal() )
      #obsZadok89.append(i[1].Zadok89.toordinal() )
    e=self.run_sites(sit_names,params)
    o=self.obs()
    #todo calcuar rmse
    pred=numpy.array([ e[k].loc[1,['Zadok65','Zadok89']] for k in e.keys()]  ) 
    targ=numpy.array([ o[k].loc[1,['Zadok65','Zadok89']] for k in e.keys()]  ) 
    sse=((pred - targ)**2).sum() 
    return sse
  def plot_optim(self,params=[]):
    #if type(params)!=list:params=list(params)
    sit_names=[]
    cal_sites=self.d[numpy.logical_and(numpy.logical_not(pandas.isna(self.d.Zadok65)),numpy.logical_not(pandas.isna(self.d.Zadok89)))]
    for i in cal_sites.iterrows():
      year=i[0][0]
      site=i[0][1]
      toss=i[0][2]
      sit_names.append(str(site)+'-'+str(year)+'-'+str(toss))
    e=self.run_sites(sit_names,params)
    o=self.obs()
    #todo calcuar rmse
    pred=numpy.array([ e[k].loc[1,['Zadok65','Zadok89']] for k in e.keys()]  ) 
    targ=numpy.array([ o[k].loc[1,['Zadok65','Zadok89']] for k in e.keys()]  ) 
    pylab.plot(pred,targ) 
    return 
  def gen_results(self,params=[46.2709270, 29.9722860,  0.0000683]):
    sit_names=[]
    cal_sites=self.d
    buf='Year\tSite\tsowDay\tsimulatedZadok10_dd/mm/yyyy\tsimulatedZadok30_dd/mm/yyyy\tsimulatedZadok65_dd/mm/yyyy\tsimulatedZadok90_dd/mm/yyyy\n'
    for i in cal_sites.iterrows():
      year=i[0][0]
      site=i[0][1]
      toss=i[0][2]
      sit_names.append(str(site)+'-'+str(year)+'-'+str(toss))
    e=self.run_sites(sit_names,params)
    for i in e:
      site,year,__=i.split('-')
      buf+=str(year)+'\t'
      buf+=str(site)+'\t'
      buf+=str(e[i].loc[1].Date)+'\t'
      buf+='NA'+'\t'
      buf+='NA'+'\t'
      buf+=str(datetime.datetime.fromordinal(e[i].loc[1].Zadok65).strftime("%d/%m/%Y") )+'\t'
      buf+=str(datetime.datetime.fromordinal(e[i].loc[1].Zadok89).strftime("%d/%m/%Y") )+'\n'
    f=open('cal2CSIRO_results numerical PG_AndresBerger.txt','w')
    f.write(buf)
    f.close()
    return 



#////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////  
from scipy.optimize import basinhopping

class MyTakeStep(object):
 def __init__(self,bounds=None): #stepsize=0.5
   #self.stepsize = stepsize
   self.bounds=bounds
 def __call__(self, x):
   #if self.bounds=None: self.bounds=[(30.,60.),(30.,45.),(-0.03,-0.003)]#,(0.003,0.033)]#,(10,50)]
   x=numpy.array([(numpy.random.uniform()*(i[1]-i[0]))+i[0] for i in self.bounds])
   return x

def opt():
  
  #MTDV=32      #Minimum thermal days for vegetative growth phase                                     d
  #MTDR=31.     #18,14#Minimum thermal days for reproductive (seed fill)phase                               d
  #PSEN=-0.06   #Photoperiod sensitivity of phenological development                                  h-1
  #VSEN= 0.003  #Vernalization sensitivity   (0.033 WW ,0.003 SWhigh latitude , 0 SW)                     d
  #VDSA=50
  #initial_guess=[40.,31.,0.03,0.0003,50]
  bounds=[(20.,80.),(20.,60.),(0.,0.3),(0.,0.1),(0.,60.)]
  s=sim()
  
  ##minimizer_kwargs = {"method": 'L-BFGS-B','bounds':bounds,'jac':'2-point'}#'options':{'finite_diff_rel_step':'2-point'}}#{'maxiter':1000}}
  #minimizer_kwargs = {"method": 'Nelder-Mead'}#L-BFGS-B TNC
  #mytakestep=MyTakeStep(bounds)
  #ret = scipy.optimize.basinhopping(s.run_optim, initial_guess,  minimizer_kwargs=minimizer_kwargs, niter=50,disp=True)#,stepsize=5)take_step=mytakestep,
  
  ret =scipy.optimize.shgo(s.run_optim, bounds,  sampling_method='sobol', n=1000,minimizer_kwargs={'method': 'Nelder-Mead'})
  #array([4.6406250e+01, 3.1484375e+01, 9.9218750e-03, 2.5781250e-03])
  #      [4.6406250e+01, 3.1484375e+01, 9.9218750e-03, 2.5781250e-03, 1.8125000e+01]
  #ret= scipy.optimize.dual_annealing(s.run_optim, bounds, x0=initial_guess)
  return ret

#croptimizr 46.2709270, 29.9722860,  0.0000683


    #Nedler 500 #38.6492883, 28.3935396, 8.54240276e-03, 7.37014968e-03, 5.68360096e+01
    #BFS 500 #3.84570312e+01, 3.72265625e+01, 7.91015625e-03, 6.54296875e-03, 5.51367188e+01
    #nedler 1000 
         #fun: 2466.0
    #funl: array([2466.])
 #message: 'Optimization terminated successfully.'
    #nfev: 2131
     #nit: 2
   #nlfev: 131
   #nlhev: 0
   #nljev: 0
 #success: True
       #x: array([4.75366505e+01, 2.88608562e+01, 1.04767787e-02, 4.15774740e-03,
       #2.58304351e+01])
      #xl: array([[4.75366505e+01, 2.88608562e+01, 1.04767787e-02, 4.15774740e-03,
        #2.58304351e+01]])


    #BFS 1000 4.75366505e+01, 2.88608562e+01, 1.04767787e-02, 4.15774740e-03,   2.58304351e+01
    
    #nedler 1000 only 3 params
         #fun: 2466.0
    #funl: array([2466.])
 #message: 'Optimization terminated successfully.'
    #nfev: 2131
     #nit: 2
   #nlfev: 131
   #nlhev: 0
   #nljev: 0
 #success: True
       #x: array([4.75366505e+01, 2.88608562e+01, 1.04767787e-02, 4.15774740e-03,
       #2.58304351e+01])
      #xl: array([[4.75366505e+01, 2.88608562e+01, 1.04767787e-02, 4.15774740e-03,
        #2.58304351e+01]])
