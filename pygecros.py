"""
Modelo trigo Jun 2014
"""
import math,sys
from pycrop import gecros_utils
#import gecros_for_phs_modiff as gecros_for
#import gecros_for_t as gecros_for
import numpy
import pylab#,rpy
#from support_data import *
#from mx import DateTime
import pandas, datetime

pylab.ion()

def print_exception(type=None, value=None, tb=None, limit=None) :
    if type is None:
      type, value, tb=sys.exc_info()
    import traceback
    output = ['Traceback (innermost last) : ' ]
    list = traceback.format_tb(tb, limit) + traceback.format_exception_only(type,value)
    for item in list :
      print (item)

"""***********************************************************************
*                               GECROS                                *
*     Genotype-by-Environment interaction on Crop growth Simulator    *
*          (as linked with a simple soil simulation model)            *
*                                                                     *
*                          (FST-version 1.0)                          *
*                         Author: Xinyou YIN                          *
*                      Crop and Weed Ecology group                    *
*                Wageningen University & Research Centre              *
*               PO Box 430, 6700 AK Wageningen, Netherlands           * 
***********************************************************************
*               Pythonized by aberger@inia.org.uy 2014                *
*                        INIA La Estanzuela                           *
*                            Colonia, UY                              *
***********************************************************************"""

def pheno(DS,SLP,DDLP,SPSP,EPSP,SVSP,EVSP,PSEN,VSEN,VDSA,CVER,MTDV,MTDR,TDU):
    """*----------------------------------------------------------------------*
    *  SUBROUTINE PHENO                                                       *
    *  Purpose: This subroutine calculates phenological development rate.     *
    *                                                                         *
    *  FORMAL PARAMETERS:  (I=input,O=output,C=control,IN=init,T=time)        *
    *                                                                         *
    *  name   type meaning                                       units  class *
    *  ----   ---- -------                                       -----  ----- *
    *  DS      R4  Development stage                               -       I  *
    *  SLP     R4  Crop type(1. for short-day,-1. for long-day)    -       I  *
    *  DDLP    R4  Daylength for photoperiodism                    h       I  *
    *  SPSP    R4  DS for start of photoperiod-sensitive phase     -       I  *
    *  EPSP    R4  DS for end of photoperiod-sensitive phase       -       I  *
    *  SVSP    R4  DS for start of vernalization-sensitive phase   -       I  *
    *  EVSP    R4  DS for end of vernalization-sensitive phase     -       I  *
    *  PSEN    R4  Photoperiod sensitivity (+ for SD, - for LD)    h-1     I  *
    *  VSEN    R4  Vernalization sensitivity                       d       I  *
    *  VDSA    R4  Vernalization saturation                        d       I  *
    *  CVER    R4  Cumulative vernalization days                   d       I  *
    *  MTDV    R4  Minimum thermal days for vegetative phase       d       I  *
    *  MTDR    R4  Minimum thermal days for reproductive phase     d       I  *
    *  TDU     R4  Daily thermal-day unit                          -       I  *
    *  DVR     R4  Development rate                                d-1     O  *
    *-------------------------------------------------------------------------*"""

    #*---determining if it is for short-day or long-day crop
    if SLP<0.:
        MOP = 18.     #minimum optimum photoperiod for long-day crop
        DLP = min(MOP,DDLP)
    else:
        MOP = 11.     #maximum optimum photoperiod for short-day crop
        DLP = max(MOP,DDLP)
    #*---effect of photoperiod on development rate
    if DS<SPSP or DS>EPSP:
        EFP = 1.
    else:
        EFP = max(0., 1.-PSEN*(DLP-MOP))
    #*---vernalization
    if CVER>=VDSA or DS<SVSP or DS>EVSP: #use the same limits as for photoperiod for simplicity
        VER=1.
    else:
        VER = max(0., 1.-VSEN*(VDSA-CVER))
    
    #*---development rate of vegetative and reproductive phases
    if DS>0. and DS<1.:
        DVR   = 1./MTDV*TDU*EFP*VER
    else:
        DVR   = 1./MTDR*TDU
    return DVR
  
def vern(TMAX,TMIN,TAVSS,CVER):
    """*-------------------------------------------------------------------------*
    *  SUBROUTINE PHENO                                                          *
    *  Purpose: This subroutine calculates the daily amount of vernalization day *
    *                                                                            *
    *  FORMAL PARAMETERS:  (I=input,O=output,C=control,IN=init,T=time)           *
    *                                                                            *
    *  name   type meaning                                    units  class       *
    *  ----   ---- -------                                    -----  -----       *
    *  TMAX    R4  Daily maximum temperature                    oC      I        *
    *  TMIN    R4  Daily minimum temperature                    oC      I        *
    *  TAVSS   R4  Soil surface temperature                     oC      I        *
    *  CVER    R4  Cumulative vernalization days                d       I        *
    *  VDU     R4  Daily amount of vernalization accumulated    d       I        *
    *----------------------------------------------------------------------------*"""
    #*---mean daily temperature
    TMP  = TAVSS#(TMAX + TMIN)/2. #should use crown temperature instead! beaware of snow cover effects not currently considered!!
    #vernalization function (Soltani & Sinclair)
    VDU=0.
    if CVER<10. and TMAX>30.:#devernalization
      VDU=-1.*(min(CVER,0.5*(TMAX-30.))) 
    else:                    #vernalization
      tbv=-1.;tp1v=0.;tp2v=8.;tcv=12.
      if   TMP<=tbv : VDU=0.
      elif TMP< tp1v: VDU= (TMP-tbv)/(tp1v-tbv)
      elif TMP<=tp2v: VDU=1.
      elif TMP< tcv : VDU=(tcv-TMP)/(tcv-tp2v)
      elif TMP>=tcv : VDU=0.
    return VDU  
  
  
  
def run(
    #********* SELF-DEFINED WATER AND NITROGEN SUPPLY TO THE CROP **********
    # User-defined daily water and (NH4+ and NO3-)nitrogen availability
    # WSWI or NSWI = 1. for using simulated soil water or nitrogen supply;
    # WSWI or NSWI =-1. for using user self-defined soil water supply (i.e.
    # WINPUT) or nitrogen supply (NINPA & NINPN). In this case, crop model
    # is de-coupled from the example soil model; in other words, simulated
    # soil water and N availabilities are no longer affecting crop growth.
    WSWI = 1., NSWI = 1.,
    WINPUT = 7.  ,#User-defined water supply to crop                                                    mm d- 1
    NINPA  = 0.0  ,#User-defined amonium-N supply to crop                                                g N m-2 d- 1
    NINPN  = 0.65 ,#User-defined nitrate-N supply to crop                                                gNm-2 d- 1
    
    #***************************** MODEL INPUTS ****************************
    #*** Crop parameters for pea (Pisum sativum L.)        
    # LEGUME = 1. for leguminous crops;   = -1. for non-leguminous crops.
    # C3C4   = 1. for C3 crops;           = -1. for C4 crops.
    # DETER  = 1. for determinate crops;  = -1. for indeterminate crops.
    # SLP    = 1. for short-day crops;    = -1. for long-day crops.
    # LODGE  = 1. for cases with lodging; = -1. for cases without lodging."""
    LEGUME = -1., C3C4 = 1., DETER = 1., SLP = -1., LODGE = -1.,   
  
    EG=1.        ,# Efficiency of germination (epsilon g)                                                gg-1
    CFV=0.48     ,# Carbon fraction in the vegetative organs  (Fcv)                                      gCg-1
    YGV=0.81     ,# Growth efficiency for vegetative organs (i.e leaves, stem roots) (Ygv)               gCg-1C
    FFAT=0.02    ,# Fraction of fat in the storage organs (Flip)                                         gfat g-1
    FLIG=0.06    ,# Fraction of lignin in the storage organs (Flig)                                      glig g-1
    FOAC=0.02    ,# Fraction of organic acids in the storage organs (Foac)                               goac g-1
    FMIN=0.02    ,# Fraction of minerals in the storage organs (Fmin)                                    gmin g-1
    LNCI=0.055   ,# 0.05Initial value of N concentration in living leaves (ncri0)                            gNg-1
    TBD=0.       ,# Base temperature for phenology (Tb)                                                  C
    TOD=25.0     ,# Optimum temperature for phenology (To)                                               C
    TCD=37.      ,# Ceiling temperature for phenology (Tc)                                               C
    TSEN=1.5     ,#Curvature for temperature response (ct)   (value range 0.5-3.0)                      --
    SPSP=0.2        ,#Development stage (DS) for start of fotoperiod sensitive phase (v1) (0-0.5)          --
    EPSP=0.7        ,#Development stage (DS) for end of photoperiod sensitive phase  (v2) (0.5-0.8)        --
    SVSP=0.        ,#Development stage (DS) for start of vernalization sensitive phase (v1) (0-0.5)          --
    EVSP=0.5        ,#Development stage (DS) for end of vernalization sensitive phase  (v2) (0.5-0.8)        --    
    INSP=-2.     ,# Inclination of sun angle for computing daylenght for photoperiodism(alpha)           degrees
    LWIDTH=0.015  ,# Leaf width (w)                                                                       m
    RDMX=100    ,# Maximum rooting depth (Dmax)                                                         cm
    KWEXT=0.12,#0.096     , #K de Daradanelli FCR 87 (2004) 59-71, tasa de extraccion de agua maxima de suelo totalmente explorado    
    CDMHT=300,#180.   ,#460. Stem dry weight per unit of plant height (rho)                                       gm-2m-1
    PMEH=0.8        ,#0.8#Fraction of sigmoid curve inflection in entire plant height growth period (0.6-0.9)  --
    ESDI=1.1        ,#Development stage for end of seed number determining period for indeterminate crops  (1.1-1.45) --
    PMES=0.6        ,#Fraction of sigmoid curve inflection in entire seed growth period (0.4-0.7)          --
    CCFIX=6.     ,# Carbon cost of symbiotic nitrogen fixation                                           g C g-1 N
    NUPTX=.45      ,#.45 Maximum crop nitrogen uptake  (0.4-0.8)                                              gN m-2 d-1
    SLA0=0.028 ,#changed from 0.028# Specific leaf area constant                                                          m2 leaf g-1
    SLNMIN=0.35  ,#0.35#changed from 0.35# Minimum or base SLN for photosynthesis                                               gNm-2
    RNCMIN=0.005 ,# Minimum N concentration in the roots                                                 gN g-1
    STEMNC=0.01 ,#0.01 Nitrogen concentration in the stem                                                   gNg-1
    WRB=0.25        ,#chg from 0.25 AGB #Critical root weight density                           g m-2 cm-1 depth
    EAJMAX=48270.,#48270.# Energy of activation for JMAX                                                      J mol-1
    XVN=45.         ,#Slope of linear relationship between VCMX and leaf nitrogen                          umol CO 2 s-1 g-1 N
    XJN=100.        ,#Slope of linear relationship between JMAX and leaf                                   umol e- S-1 g-1 N
    THETA=0.8       ,#Convexity for light response of electron transport (J2) in photosynthesis            --
    #*** Genotype-specific parameters 
    SEEDW=0.035 ,#0.0280Seed weight                                                                          g seed-1
    SEEDWF=0.18, #Proportion of seed usable for initial plant size
    SEEDNC=0.0272,#0.0272,#14.5%PCbs 0.024 Standard seed (storage organ) N concenh'ation                    g N g-1
    BLD=45.      ,#0.40,# Leaf angle (from horizontal)                                                         degree
    HTMX=0.8    ,# Maximum plant height                                                                 m
    MTDV=32,   #Minimum thermal days for vegetative growth phase                                     d
    MTDR=31.,    #18,14#Minimum thermal days for reproductive (seed fill)phase                               d
    PSEN=-0.06,       #Photoperiod sensitivity of phenological development                                  h-1
    VSEN= 0.003,       #Vernalization sensitivity   (0.033 WW ,0.003 SWhigh latitude , 0 SW)                     d
    VDSA=50      ,#winterW USA #Vernalization saturation day units                                       d
    #*** Soil parameters
    PNLS=0.5      ,# Fraction of dead leafN incorporated into soil litter N   
    CLAY=23.4, WCMIN=0.2, WCFC=0.45, WCMAX=0.54,
    DRPM=1.44       ,#Ratio DPM/RPM of added plant material                                                --
    DPMR0=10.       ,#Standard value for DPMR                                                              yr-1
    RPMR0=0.3       ,#Standard value for RPMR                                                              yr-1
    BIOR= .66      ,#0.66Decomposition rate constant for BIO                                                  yr-1
    HUMR=0.02       ,#0.02Decomposition rate constant for HUM                                                  yr-1
    TOC= 20000,#1300000*0.014=18000#7193   ,#7193.#Total organic C in the soil                                                          g C m-2
    BHC= 3000      ,#4400.#Initial value for (BIO + HUM)                                                        g C m-2
    FBIOC=0.03      ,#Fraction of BIOI in initial total soil organic carbon (TOC)                          --
    RN=1.5, RA=1.5    ,#Residual N in the soil (not plant uptakeable)                                        g N m-2
    RSS=100.        ,#Soil resistance, equivalent to leaf stomatal resistance                              s m-1
    SD1=15.          ,#10#initial roothing depth                                                               cm 
    TCT=1.          ,#Time constant for soil temperature dynamics                                          d
    TCP=1.          ,#Time constant for some soil dynamic processes (-1 d)                                 d
    MULTF=1.05        ,#Multiplication factor for initial soil water status (percent of FC)                  --
    #*** Sensitivity-analysis options
    CO2A=385.     ,# Ambient CO2 concentration                                                           umol mol-1 
    COEFR=1.      ,# Factor for change in radiation, for sensitivityanalysis                              --
    COEFV=1.      ,# Factor for change in vapour pressure, for sensitivity analysis                       --
    COEFT=0.      ,# Increment of a temperature, for sensitivity analysis                                 --
    FCRSH=0.5,#5     ,# 0.5 Initial fraction of carbon in the shoot                                              --
    FNRSH=0.63,#0.63#6    ,# 0.63  Initial fraction ofN in the shoot                                                    --
    PNPRE=0.8     ,# 0.7 Propotiion of seed N that comes from non-structural N in vegetative organs(0.6-0.95) --
                    # accumulated before end of seed-number determining period                             --
    CB=0.75       ,#Factor for initial N concentration of seed fill                                      --
    CX=1.         ,#Factor for final N concentration of seed fill                                        --
    TM=1.5        ,#DS when transition from CB to CX is fastest                                          --
    #*** Management actions                
    NPL    = 250. ,#Plant density                                                                        plants m-2
    STTIME ='2011-6-20'   ,#start time                                                                           
    FINTIM ='2020-12-10'   ,#time to finish run                                                                   
    IRRIG  = [[5,0],[15,0],[25,0],[35,0],[45,0]] , #Irrigation schedule                                  mm d-1
    FERTNA = [[5,0],[15,0],[25,0],[35,0],[45,0]] , #Fertilizer amonium schedule                          g N m-2d-1
    FERTNN = [[5,0],[15,0],[25,0],[35,0],[45,0]] , #Fertilizer nitrate schedule                          g N m-2d-1
    NAI    = 6.,#1.5   ,#Initial value for NA (amonium N in the soil available for uptake)                    g N m-2
    NNI    = 14.,#3.   ,#Initial value for NN (nitrate N in the soil available for uptake)                    g N m-2
    WEATHER=None,
    LAT=34.3
    ):
    
    def afgen(F,X):
        for i in range(len(F)-1):
            if F[i+1][0]>X:break
        if F[i][0]>X:return X
        return F[i][1]+((X-F[i][0])*(F[i+1][1]-F[i][1]))/(F[i+1][0]-F[i][0])
    #************************* RUN CONTROL *********************************
    DELT=1.#time step in days
    PRDEL=1.#time step for output in days  
    
    YEAR   =STTIME.year
    STTIME =STTIME.toordinal()   #start time                                                                           
    FINTIM =FINTIM.toordinal()   #time to finish run                                                                   
        
    ## create time series variables
    d={}
    
    #************************ INITIAL CONDITIONS ***************************
    ZERO   = 0.

    #*** Initial conditions for crop model
    FPRO  = 6.25*SEEDNC
    FCAR  = 1.-FPRO-FFAT-FLIG-FOAC-FMIN
    CFO   = 0.444*FCAR+0.531*FPRO+0.774*FFAT+0.667*FLIG+0.368*FOAC
    YGO   = CFO/(1.275*FCAR+1.887*FPRO+3.189*FFAT+2.231*FLIG+ \
            0.954*FOAC)*30./12.

    CLVI   = NPL * SEEDWF * SEEDW * CFO * EG * FCRSH
    CRTI   = NPL * SEEDWF * SEEDW * CFO * EG * (1.-FCRSH)
    
    NLVI   = LNCI* CLVI/CFV
    NRTI   = NPL * SEEDWF * SEEDW * EG * LNCI * FCRSH/FNRSH - NLVI

    LNCMIN = SLA0*SLNMIN
    LAII   = CLVI/CFV*SLA0
    SLNBI  = NLVI/LAII
    
    RDMX1=RDMX+1. #in original model fixed to RDMX1=150.
    RDI    = max(2., SD1)
    HTI    = HTMX/1000.

    #*** Initial conditions for the example soil model**********************

    TSOILI = 15.
    WCI    = WCFC * MULTF
    WULI   = 10.*(WCI-WCMIN)*RDI
    WLLI   = 10.*(WCI-WCMIN)*(RDMX1-RDI)

    DPMI   = ZERO
    RPMI   = TOC - BHC - DPMI
    BIOI   = FBIOC * TOC
    HUMI   = BHC - BIOI

    DPNI   = 1./ 40.*DPMI
    RPNI   = 1./100.*RPMI
    
    NAULI  = (1.-math.exp(-0.065*RDI))*(NAI) +     RDI/RDMX1 *RA #(NAI-RA) was before (NAI)
    NALLI  =     math.exp(-0.065*RDI) *(NAI) + (1.-RDI/RDMX1)*RA 
    NNULI  = (1.-math.exp(-0.065*RDI))*(NNI) +     RDI/RDMX1 *RN 
    NNLLI  =     math.exp(-0.065*RDI) *(NNI) + (1.-RDI/RDMX1)*RN 
    
    #*** Initialize stock variables ****************************************
    DS=ZERO    #Development stage                                                                        --
    CTDU=ZERO  #Cummulative thermal day unit                                                             d
    CVER=ZERO  #Cumulative vernalization days                                                            d
    
    CLV=CLVI   #Amount of carbon in the living leaves                                                    gCm-2
    CLVD=ZERO  #Amount of carbon in the dead leaves                                                      gCm-2
    CSST=ZERO  #Amount of carbon in structural stems                                                     gCm-2 
    CSO=ZERO   #Amount of carbon in storage organs                                                       gCm-2
    CSRT=CRTI  #Amount of carbon in structural roots                                                     gCm-2
    CRTD=ZERO  #Amount of carbon in dead roots                                                           gCm-2
    CLVDS=ZERO #Amount of carbon in the dead leaves that had become litters in the soil                  gCm-2
    NRT=NRTI   #Nitrogen in living roots                                                                 g N m- 2
    NST=ZERO   #Nitrogen content in stems (including structural stem and reserves)                       gNm-2 
    NSO=ZERO   #Nitrogen content in storage organs                                                       gNm-2
    TNLV=NLVI  #Total leaf nitrogen (including N in senesced leaves)                                     g N m-2
    NLVD=ZERO  #Nitrogen in dead leaves                                                                  g N m- 2
    NRTD=ZERO  #Nitrogen in dead roots                                                                   g N m-2

    NLV=NLVI   #Nitrogen in living leaves                                                                gNm-2
    
    CRVS=ZERO  #Amount of carbon in stems reserves                                                       gCm-2
    CRVR=ZERO  #Amount of carbon in root reserves                                                        gCm-2
    NREOE=ZERO #NRES accumulated till the end of seed-number determining period                          g N m-2
                #NRES == Estimated vegetative-organ N remobilizable for seed growth
    NREOF=ZERO #NRES accumulated till the the moment at which seed fill starts                           g N m-2
    DCDSR=ZERO #Shortfall of carbon demand for seed fill in previous time steps                          g C m-2
    DCDTR=ZERO #Shortfall of carbon demand for structural stems in previous time steps                   g C m- 2
    SLNB=SLNBI #SLN in bottom leaves of canopy                                                           gNm-2 leaf
    LAIC=LAII  #Carbon determined LAI  m2leafm-2                                                         m2leafm-2
    LAI =LAII
    RMUL=ZERO  #Total ofRMUN+ RMUA+RMUS+RMLD                                                             g CO 2 m-2 d-1
    NDEMP=ZERO #NDEM of the previous time step                                                           g N m-2 d-1
    NSUPP=ZERO #NSUP of the previous time step                                                           g N m-2 d-1
    NFIXT=ZERO #Total symbiotically fixed nitrogen during growth                                         g N m-2
    NFIXR=ZERO #Reserve pool of symbiotically fixed nitrogen                                             g N m-2
    DCDTP=ZERO #Carbon demand for structural stem growth at the previous time step                       g C m-2 d-1
    HT=HTI     #Plant height                                                                             m
    RD=RDI     #Rooting depth to the soil                                                                cm
    TDAPAR=ZERO#Total PAR absorbed by the canopy during growth                                           J m-2
    TPCAN=ZERO #Cumulative canopy photosynthesis over growth period                                      g CO2 m-2
    TRESP=ZERO #Total crop respiratory cost during growth                                                g CO2 m-2
    TTCAN=ZERO #Cumulative canopy transpiration                                                          mm
    TNUPT=ZERO #Total crop nitrogen uptake during growth                                                 g N m-2
    LITNT=ZERO #Total LITN (Litter nitrogen entering the soil) during growth                             g N m-2

    TSOIL=TSOILI#Soil temperature                                                                        degC
    TAVSS=TSOILI # Soil surface temperature                                                              degC
    WUL=WULI   #Water content in the upper soil layer                                                    mm
    WLL=WLLI   #Water content in the lower soil layer                                                    mm
    CETA=ZERO  #Cummulative ETA                                                                          mm
    WUG=ZERO   #Cummulative water flow to groundwater                                                    mm
    DPM=DPMI   #Decomposable plant material                                                              g C m-2
    RPM=RPMI   #Resistant plant material (difficult to decompose)                                        g C m-2
    BIO=BIOI   #Microbial biomass in the soil                                                            g C m-2
    HUM=HUMI   #Humified organic matter in the soil                                                      g C m-2
    DPN=DPNI   #Organic N in DPM                                                                         g N m-2
    RPN=RPNI   #Organic N in RPM                                                                         g N m-2
    
    TMDN=ZERO  #Cummulative mineralized N in the whole soil layer                                        g N m-2
    TNLEA=ZERO #Total nitrate N leached to groundwater                                                   g N m-2
    NAUL=NAULI #Al11111onium-N in the upper soil layer                                                   g N m-2
    NALL=NALLI #Ammonium-N in the lower soil layer                                                       g N m-2
    NNUL=NNULI #Nitrate-N in the upper layer                                                             g N m-2
    NNLL=NNLLI #Nitrate-N in the lower layer                                                             g N m-2
    SFERNA=ZERO #amonium N feliilizer susceptible to volatilization                                      g N m-2

    #load WEATHER data
    #if WEATHER==None:
      ##WEATHER=numpy.fromfile(file='weather_data.txt',sep=' ')
      #WEATHER=numpy.fromfile(file='wheather_data_LE1965-2014.txt',sep=' ')
      ##WEATHER=numpy.fromfile(file='Libro2.txt',sep='\t')
      #WEATHER=WEATHER.reshape((numpy.shape(WEATHER)[0]/9,9)).T
    #WEATHER_POINTER=0
    #</INITIAL>
    anthesis=False
    TIME=STTIME
    print ('emergencia ',datetime.datetime.fromordinal(TIME).isoformat() )
    print (DS > 2. ,TIME>FINTIM)
    fail=False
    while not (DS > 2. or TIME>FINTIM or fail):
        #print ('*************** DAY %s ********************'%(TIME))
        #********************* ENVIRONMENTAL DATA ******************************
        
        try:
          #WEATHER CNTR='NLD'; ISTN=1; WTRDIR ='D:\WEATHER\'; IYEAR=2003
          #*     RDD    Daily global radiation in     J/m2/d
          #*     TMMN   Daily minimum temperature in  degree C
          #*     TMMX   Daily maximum temperature in  degree C
          #*     VP     Vapour pressure in            kPa
          #*     WN     Wind speed in                 m/s
          #*     RAIN   Precipitation in              mm
          #*     LAT    Latitude of the side          degree
          #*     DOY    Day of year                   d
          
          #print (TIME,STTIME)
          if 0:
            # WEATHER file rows:
            #    0 station number
            #    1 year
            #    2 day
            #    3 irradiation (kJ m-2 d-1)
            #    4 minimum temperature (degrees Celsius)
            #    5 maximum temperature (degrees Celsius)
            #    6 vapour pressure (kPa)
            #    7 mean wind speed (m s-1)
            #    8 precipitation (mm d-1)
            t=DateTime.DateTimeFromJDN(TIME)
            for i in range(WEATHER_POINTER,len(WEATHER[2])):
                if int(WEATHER[1][i])==t.year and int(WEATHER[2][i])==t.day_of_year:
                    WEATHER_POINTER=i
                    RDD=WEATHER[3][i]*1000.
                    TMMN=WEATHER[4][i]
                    TMMX=WEATHER[5][i]
                    VP=WEATHER[6][i]
                    WN=WEATHER[7][i]
                    RAIN=WEATHER[8][i]
                    #LAT=-34.5
                    DOY=t.day_of_year
                    #print (WEATHER[1][i],YEAR,weather[2][i],TIME#'GOT WEATHER !!')
                    break
          else:
            #AGMIP
            #'SRAD', u'TMAX', u'TMIN', u'RAIN', u'WIND', u'DEWP', u'VPRS', u'RHUM'
            #t=int(str(int(DateTime.DateTimeFromJDN(TIME).year))+'{number:0{width}d}'.format(width=3, number=int(DateTime.DateTimeFromJDN(TIME).day_of_year)) )
            t= pandas.Timestamp(datetime.datetime.fromordinal(TIME ))
            #print (TIME)
            #print (t)
            #print (WEATHER[2])
            RDD= WEATHER[2]['SRAD'][t:t].values[0]*1e6    
            TMMN=WEATHER[2]['TMIN'][t:t].values[0]
            TMMX=WEATHER[2]['TMAX'][t:t].values[0]    
            VP=  WEATHER[2]['VPRS'][t:t].values[0]*10.    
            WN=  WEATHER[2]['WIND'][t:t].values[0]    
            RAIN=WEATHER[2]['RAIN'][t:t].values[0]    
            DOY=datetime.datetime.fromordinal(TIME ).timetuple().tm_yday#DateTime.DateTimeFromJDN(TIME).day_of_year
            #LAT=45.77
            #print (RDD,TMMN,TMMX,VP,WN,RAIN,DOY)
        except:
          print_exception()
          fail=True
        DFS    = int(TIME - STTIME) #days from sstime as int
        
        if not anthesis and DS>1.:
          print ('antesis',DFS,datetime.datetime.fromordinal(TIME).isoformat())
          anthesis=True
          
        TMAX   = TMMX + gecros_utils.insw (DS-0., 0., COEFT)
        TMIN   = TMMN + gecros_utils.insw (DS-0., 0., COEFT)
        DAVTMP = 0.29*TMIN + 0.71*TMAX
        NAVTMP = 0.71*TMIN + 0.29*TMAX
    
        DDTR   = RDD * gecros_utils.insw(DS-0., 0., COEFR)
        DVP    = VP  * gecros_utils.insw(DS-0., 0., COEFV)
        WNM    = max(0.1, WN) 
        
        #********************** THE EXAMPLE SOIL MODEL *************************
        #*** Soil water availability and water balance
        IRRI=0;FERNA=0;FERNN=0
        for i in IRRIG:
            if i[0]==TIME: IRRI=i[1]
        for i in FERTNA:
            if i[0]==TIME: 
                FERNA=i[1] 

        for i in FERTNN:
            if i[0]==TIME: FERNN=i[1]  
        RFIR   = RAIN + IRRI
        
        #*** Soil organic nitrogen
        CNDRPM = (DPM+RPM)/gecros_utils.notnul(DPN+RPN)
        
        #*** Soil organic carbon
        CBH    = 1.67*(1.85+1.60*math.exp(-0.0786*CLAY))
        FT     = 47.9/(1.  +math.exp(106./max(0.2,TSOIL+18.3)))#0.2 min to avoid float overflow
        FM     = gecros_utils.limit(0.2, 1.0, 0.2+0.8*(WUL+WLL)/10./RDMX1/(WCFC-WCMIN))   
        
        DPMRC  = gecros_utils.insw(NNUL+NAUL+NNLL+NALL-RA-RN, 0., DPMR0)
        DPMR   = gecros_utils.insw(1./gecros_utils.notnul(CNDRPM)-1./8.5/(1.+CBH), DPMRC, DPMR0)
        RPMRC  = gecros_utils.insw(NNUL+NAUL+NNLL+NALL-RA-RN, 0., RPMR0) 
        RPMR   = gecros_utils.insw(1./gecros_utils.notnul(CNDRPM)-1./8.5/(1.+CBH), RPMRC, RPMR0)
  
                          
        DECDPM = DPM*(1.-math.exp(-FT*FM*DPMR/365.))/TCP
        DECRPM = RPM*(1.-math.exp(-FT*FM*RPMR/365.))/TCP
        DECBIO = BIO*(1.-math.exp(-FT*FM*BIOR/365.))/TCP
        DECHUM = HUM*(1.-math.exp(-FT*FM*HUMR/365.))/TCP       
              
        

        RESCO2 = CBH /(1.+CBH)*(DECDPM+DECRPM+DECBIO+DECHUM)
    
        #*** Soil organic nitrogen
        DECDPN = DPN*(1.-math.exp(-FT*FM*DPMR/365.))/TCP
        DECRPN = RPN*(1.-math.exp(-FT*FM*RPMR/365.))/TCP    
                
    
        #*** Soil mineral nitrogen
        NA     = NAUL + NALL
        NN     = NNUL + NNLL
        NMINER = NA   + NN
    
        MDN    = 1./8.5*(DECBIO+DECHUM)+ DECDPN+DECRPN - \
                    1./8.5/(1.+CBH)*(DECDPM+DECRPM+DECBIO+DECHUM)
        
        MDNUL  = (1.-math.exp(-0.065*RD))*MDN
        MDNLL  =     math.exp(-0.065*RD) *MDN
        MINAUL = gecros_utils.insw(MDN,-min((NAUL-      RD /RDMX1*RA)/TCP,-MDNUL),MDNUL)
        MINALL = gecros_utils.insw(MDN,-min((NALL-(RDMX1-RD)/RDMX1*RA)/TCP,-MDNLL),MDNLL)
        MINNUL = gecros_utils.insw(MDN,-min(NNUL/TCP,-MDNUL+MINAUL), 0.)
        MINNLL = gecros_utils.insw(MDN,-min(NNLL/TCP,-MDNLL+MINALL), 0.)
        
        FMUL   = gecros_utils.limit(0.2, 1.0, 0.2+0.8*WUL/10./      RD /(WCFC-WCMIN))
        FMLL   = gecros_utils.limit(0.2, 1.0, 0.2+0.8*WLL/10./(RDMX1-RD)/(WCFC-WCMIN))
        NITRUL = max(0.,(NAUL+MINAUL*TCP-      RD /RDMX1*RA))*  \
                    (1.-math.exp(-FT*FMUL*0.6/7.))/TCP
        NITRLL = max(0.,(NALL+MINALL*TCP-(RDMX1-RD)/RDMX1*RA))*  \
                    (1.-math.exp(-FT*FMLL*0.6/7.))/TCP
        
        DENIUL = .0005*max(0.,NNUL+MINNUL*TCP-      RD /RDMX1*RN)*  \
                    RESCO2*(1.-math.exp(-0.065*RD))
        DENILL = .0005*max(0.,NNLL+MINNLL*TCP-(RDMX1-RD)/RDMX1*RN)*  \
                    RESCO2*    math.exp(-0.065*RD)
                    
        FWS    = min(1., WUL/(RD*10.*(WCFC-WCMIN)))
        
        NSUPAS = max (0., NAUL+(MINAUL-NITRUL)*TCP-RD/RDMX1*RA)/TCP#TODO esto deberia ser ponderado por la exploracion radicular, no se explora toda la superficie, haerlo en funcion de LAI como en campbell&diaz
        NSUPNS = max (0., NNUL+(MINNUL-DENIUL)*TCP-RD/RDMX1*RN)/TCP*FWS
        
        #Variables pased to plant model
        
        NSUPA  = gecros_utils.insw(NSWI, NINPA, NSUPAS)#*min(1.,max(0.05,1-math.exp(-0.4*LAI)))
        NSUPN  = gecros_utils.insw(NSWI, NINPN, NSUPNS)#*min(1.,max(0.05,1-math.exp(-0.4*LAI)))
        NSUP   = NSUPA + NSUPN
        DWSUP  = gecros_utils.insw(WSWI, WINPUT, max(0.1,KWEXT*WUL/TCP+0.1))
        WCUL   = (WUL+WCMIN*10.*RD)/10./RD
        

        #********************** THE PART OF CROP DYNAMICS **********************

        #*** Photoperiod, solar constant and daily extraterrestrial radiation
        
        SC,SINLD,COSLD,DAYL,DDLP,DSINBE=gecros_utils.astro(DOY,LAT,INSP)

        DLAI   = (CLVD-CLVDS)/CFV*SLA0
        TLAI   = LAIC + CLVD /CFV*SLA0
        TDLAI  = (CLVD+CLVDS)/CFV*SLA0
        
        #*** Extinction coefficient of nitrogen and wind
        KL=gecros_utils.kdiff(TLAI,BLD*3.141592654/180.,0.2)
    
        KLN    = KL*(TNLV-SLNMIN*TLAI)
        NBK    = SLNMIN*(1.-math.exp(-KL*TLAI))
        KN     = 1./TLAI*math.log((KLN+NBK)/(KLN*math.exp(-KL*TLAI)+NBK))
        KW     = KL

        #*** Leaf area development
        LAIN   = math.log(1.+KN*max(0.,NLV)/SLNMIN)/KN
        LAI    = min(LAIN, LAIC) #if LAIN>LAIC hay sensescencia luego...
            
        
        #*** Specific leaf nitrogen and its profile in the canopy
        if LAI==0:LAI=max(LAIN, LAIC)
        SLN    = NLV/LAI
        SLNT   = NLV*KN                  /(1.-math.exp(-KN*LAI))
        SLNBC  = NLV*KN*math.exp(-KN*LAI)/(1.-math.exp(-KN*LAI))
        SLNNT  = (NLV+0.001*NLV)*KN /(1.-math.exp(-KN*LAI))
        RSLNB  = (SLNBC-SLNB)/DELT
        
        VLS  = [(0.,0.),(2.5,0.)]
        LS     = gecros_utils.insw(LODGE, 0., afgen(VLS,DS))
        #FUNCTION VLS  = 0.,0., 2.5,0.
        FVPD   = gecros_utils.insw (C3C4, 0.195127, 0.116214)
        
        PPCAN,APCANS,APCANN,APCAN,PTCAN,ATCAN,PESOIL,AESOIL,\
        DIFS,DIFSU,DIFSH,DAPAR,CANT,CANMAX=gecros_utils.totpt(SC,SINLD,COSLD,DAYL,DSINBE,DDTR,TMAX,TMIN,DVP, \
                      WNM,C3C4,LAIC,TLAI,HT,LWIDTH,RD,SD1,RSS,BLD,KN,    \
                      KW,SLN,SLNT,SLNNT,SLNMIN,DWSUP,CO2A,LS,EAJMAX,     \
                      XVN,XJN,THETA,WCUL,FVPD) 
        
        #*** Developmental stage (DS) & cumulative thermal units (CTDU)
        TDU=gecros_utils.tunit(DS,TMAX,TMIN,max(0.,DIFS),DAYL,TBD,TOD,TCD,TSEN)
        VDU=vern(TMAX,TMIN,TAVSS,CVER)
        DVR=pheno(DS,SLP,DDLP,SPSP,EPSP,SVSP,EVSP,PSEN,VSEN,VDSA,CVER,MTDV,MTDR,TDU)
        #DVR=gecros_utils.pheno(DS,SLP,DDLP,SPSP,EPSP,PSEN,MTDV,MTDR,TDU)
        #DVR=pheno(DS,a,b,c,DDLP,MTDR,TMAX,TMIN,TDU)
        #*** Nitrogen partitioning between shoots and roots
        NSH    = NST + NLV + NSO
        NTOT   = NSH + NRT 
        
        #*** Biomass formation
        WLV    = CLV  / CFV
        WST    = CSST / CFV + CRVS/0.444
        WSO    = CSO  / CFO
        WRT    = CSRT / CFV + CRVR/0.444
        WSH    = WLV  + WST + WSO
        WLVD   = CLVD / CFV
        WSHH   = WSH  + (WLVD-CLVDS/CFV)#shoot includding dead minus decomposed leaves already in litter
        WSHD   = WSH  + WLVD #shoot includding dead, calculated not in original gecros
        WTOT   = WSH  + WRT
        
        HI     = WSO  / WSHD   #originaly calculated over WSHH        
        WRTD   = CRTD / CFV
              
        #*** Nitrogen concentration of biomass
        LNC    = NLV / WLV
        RNC    = NRT / WRT
        HNC    = NSH / WSH
        HNCD   = 100*NSHD/ WSHD #calculated not in original gecros
        PNC    = NTOT/ WTOT
        ONC    = gecros_utils.insw(-WSO, NSO/gecros_utils.notnul(WSO), 0.)            

        #*** Maintenance and total respiration (g CO2 m-2 d-1)
        RMRE   = max(min(44./12.*0.218*(NTOT-WSH*LNCMIN-WRT*RNCMIN),  \
                        APCAN-1.e-5-RMUL), 0.)
        RM     = max(0., min(APCAN-1.e-5,RMUL) + RMRE)
  
        #*** Nitrogen fixation (g N m-2 d-1)
        NFIXE  = max(0., APCAN-1.e-5-RM)/CCFIX*12./44.
        NFIXD  = max(0., NDEMP - NSUPP)
        NFIX   = gecros_utils.insw (LEGUME, 0., min(NFIXE, NFIXD))
        RX     = 44./12.*(CCFIX*NFIX)       
        
        #*** Current photo-assimilates (g CO2 m-2 d-1) for growth, and R/P ratio
        ASSA   = APCAN - RM - RX
        
        #*** Carbon accumulation
        CSH    = CLV + CSST+CRVS + CSO #Amount of carbon in living shoot organs (including stem reserves) gCm-2
        CRT    = CSRT+ CRVR #Amount of carbon in living roots organs (including root reserves) gCm-2
        CTOT   = CSH + CRT #Total amount of carbon in living shoots and roots gCm-2

        #*** Crop nitrogen demand and uptake (g N m-2 d-1)
        SHSA   = 12./44. * YGV*(APCAN -RM -RX)/ CSH
        RMN    = max(0., min(APCAN-1.e-5,RMUL) + max(min(44./12.*0.218*  \
                  (1.001*NTOT-WSH*LNCMIN-WRT*RNCMIN),APCAN-1.e-5-RMUL), 0.))
        SHSAN  = 12./44. * YGV*(APCANN-RMN-RX)/ CSH
        DERI   = max(0.,(SHSAN - SHSA)/(0.001*NTOT/CTOT))
        NDEMA  = CRT * (SHSA**2.)/gecros_utils.notnul(DERI)
        
        #*** Nitrogen partitioning between shoots and roots
        NCR    = gecros_utils.insw(SLNT-SLNMIN,0.,min(NUPTX,NDEMA))/(YGV* \
                (APCANS-RM-RX)*12./44.)
        
        #*** Estimation of total seed number, and 1000-seed weight
        ESD    = gecros_utils.insw(DETER, ESDI, 1.)   
        NRES   = NREOF + (NREOE-NREOF)*(ESD-1.)/gecros_utils.notnul(min(DS,ESD)-1.)
        
        TSN    = NRES/PNPRE/SEEDNC/SEEDW #original
        #TSN    = (NRES*PNPRE)/(SEEDNC*SEEDW)#corrected similar to (NRES*PNPRE)/SEEDNC/SEEDW
        TSW    = WSO/gecros_utils.notnul(TSN)*1000.
        
        #*** Carbon partitioning among organs and reserve pools
        FCSH   = 1./(1.+NCR*DERI/SHSA)
        
        #*** Carbon supply from current photo-assimilates for shoot & root growth
        DCSS   = 12./44.*    FCSH *ASSA
        DCSR   = 12./44.*(1.-FCSH)*ASSA          
        
        #*** Daily carbon flow for seed filling
        FDS=gecros_utils.betaf(DVR,1.,PMES*1.,gecros_utils.limit(1.,2.,DS)-1.)
        DCDSC,DCDS,FLWCS=gecros_utils.sinkg(DS,1.,TSN*SEEDW*CFO,YGO,FDS,DCDSR,DCSS,DELT)
        
        #*** Daily carbon flow for structural stem growth
        DCST   = DCSS - FLWCS
        IFSH   = gecros_utils.limit(0.,1.,DCST/gecros_utils.notnul(DCDTP))
        FDH=gecros_utils.betaf(DVR,(1.+ESD)/2.,PMEH*(1.+ESD)/2.,min((1.+ESD)/2.,DS))
        DCDTC,DCDT,FLWCT=gecros_utils.sinkg(DS,0.,CDMHT*HTMX*CFV,YGV,FDH*IFSH,DCDTR,DCST,DELT)
        FCSST  = gecros_utils.insw(DS-(ESD+0.2), FLWCT/DCSS, 0.)   
        
        #*** Root senescence
        KCRN   = -math.log(0.05)/6.3424/CFV/WRB/RDMX
        CSRTN  = 1./KCRN*math.log(1.+KCRN*max(0.,(NRT*CFV-CRVR*RNCMIN))/RNCMIN)
        LCRT   = max(min(CSRT-1.e-4,CSRT-min(CSRTN,CSRT)),0.)/DELT
        LWRT   = LCRT/CFV
        LNRT   = LWRT*RNCMIN  
        
        #*** Daily carbon flow for seed filling
        FCSO   = FLWCS/DCSS
        FCLV   = gecros_utils.reaand(LAIN-LAIC,ESD-DS)*(1.-FCSO-FCSST)
        
        FCRVS  = 1. - FCLV - FCSO - FCSST
        FCRVR  = gecros_utils.insw(CSRTN-CSRT, 1., 0.)
  
        #*** Leaf senescence
        LWLVM  = (LAIC-min(LAIC,LAIN))/SLA0/DELT
        LWLV   = min(WLV-1.e-5, LWLVM+gecros_utils.reanor(ESD-DS,LWLVM)*0.03*WLV) #rate of loss of leaf weight because of scenescence gm-2d-1
        LCLV   = LWLV*CFV
        LNLV   = min(LWLV,LWLVM)*LNCMIN + (LWLV-min(LWLV,LWLVM))*LNC           
        
        #*** Dynamics of carbon-reserve pool in stems and roots
        GAP    = max(0., DCDS-DCSS)
        CREMSI = min(0.94*CRVS, CRVS/gecros_utils.notnul(CRVS+CRVR)*GAP)/0.94
        CREMRI = min(0.94*CRVR, CRVR/gecros_utils.notnul(CRVS+CRVR)*GAP)/0.94
        CREMS  = gecros_utils.insw(DCDS-DCSS, 0., CREMSI)
        RCRVS  = FCRVS*DCSS - CREMS
        CREMR  = gecros_utils.insw(DCDS-DCSS, 0., CREMRI)
        RCRVR  = FCRVR*DCSR - CREMR            
        
        #*** Carbon production rate            
        RCLV   = 12./44.*ASSA*    FCSH *    FCLV  *YGV - LCLV
        RCSST  = 12./44.*ASSA*    FCSH *    FCSST *YGV
        RCSRT  = 12./44.*ASSA*(1.-FCSH)*(1.-FCRVR)*YGV - LCRT
        RCSO   = 12./44.*ASSA*FCSH*FCSO*YGO + 0.94*(CREMS+CREMR)*YGO                     
        
        #Biomass formation rates
        RWST   = RCSST/ CFV + RCRVS/0.444
        RWSO   = RCSO / CFO
        RWLV   = RCLV / CFV     #Rate of change in dry weight of living leaves gm-2d-1
        RWRT   = RCSRT/ CFV + RCRVR/0.444

        RG     = 44./12.*((1.-YGV)/YGV*(RCLV+RCSST+RCSRT+LCLV+LCRT)+  \
                            (1.-YGO)/YGO* RCSO)
        RESTOT = RM+RX+RG + 44./12.*0.06*(CREMS+CREMR)
        RRP    = RESTOT / APCAN
        
        #*** Nitrogen partitioning between shoots and roots
        NSHH   = NSH +(WLVD-CLVDS/CFV)*LNCMIN
        NSHD   = NSH + NLVD# in shoot including dead leaves and N incorporated into litter, calculated not in original gecros
        #*** Soil temperature
        TAVSS  = ((DAVTMP+DIFS)+NAVTMP)/2.
        RTSOIL = (TAVSS - TSOIL)/TCT
        
        #*** Daily and total C and N returns from crop to soil
        LVDS   = (CLVD-CLVDS)/10.*(TAVSS-TBD)/(TOD-TBD)
    
        LITC   =  LCRT + LVDS
        LITN   =  LNRT + LVDS/CFV *LNCMIN*PNLS
    
        NRETS  = LITNT+gecros_utils.insw(DS-2.,0.,NLV+NST+NRT+NFIXR+       \
                    (CLVD-CLVDS)/CFV*LNCMIN*(1.+PNLS)/2.)   
        
        #*** Soil organic nitrogen    
        RDPN   = LITN/(1.+ 40./DRPM/100.) - DECDPN
        RRPN   = LITN/(1.+100.*DRPM/40. ) - DECRPN   
        #*** Soil organic carbon
        RDPM   = LITC*DRPM/(1.+DRPM) - DECDPM
        RRPM   = LITC*1.  /(1.+DRPM) - DECRPM
        RBIO   = 0.46/(1.+CBH)*(DECDPM+DECRPM+DECBIO+DECHUM) - DECBIO
        RHUM   = 0.54/(1.+CBH)*(DECDPM+DECRPM+DECBIO+DECHUM) - DECHUM

        #*** Nitrogen partitioning between shoots and roots
        FNSH   = 1./(1.+NCR*DERI/SHSA*CSH/CRT*NRT/NSH)     

        #*** Crop nitrogen demand and uptake (g N m-2 d-1)
        HNCCR  = LNCI*math.exp(-0.4*DS)
        
        NDEMD  = gecros_utils.insw(DS-(ESD+0.8), WSH*(HNCCR-HNC)*(1.+NRT/NSH)/DELT, 0.)#Deficiency driven demand allowed in post-anthesis up to ESD+0.7
        #NDEMD  = gecros_utils.insw(DS-ESD, WSH*(HNCCR-HNC)*(1.+NRT/NSH)/DELT, 0.)#Deficiency driven demand allowed in post-anthesis up to ESD+0.7
        #NDEMD  =  WSH*(HNCCR-HNC)*(1.+NRT/NSH)/DELT #Deficiency driven demand allowed in post-anthesis
        
        NDEMAD = gecros_utils.insw(LNC-1.5*LNCI, max(NDEMA, NDEMD), 0.)
        NDEM   = gecros_utils.insw(SLNMIN-SLN+1.e-5, min(NUPTX,NDEMAD), 0.)
        NUPTA  = min(NSUPA, NSUPA/gecros_utils.notnul(NSUP)*max(0.,NDEM-NFIXR/TCP))
        NUPTN  = min(NSUPN, NSUPN/gecros_utils.notnul(NSUP)*max(0.,NDEM-NFIXR/TCP))
        NUPT   = max(0., NUPTA + NUPTN + min(NDEM, NFIXR/TCP))
        
        RNDEMP = (NDEM - NDEMP)/DELT
        RNSUPP = (NSUP - NSUPP)/DELT
        RNFIXR = NFIX - min(NDEM,NFIXR/TCP)            
        
        #*** Nitrogen accumulation rates
        #CALL 
        RNRT,RNST,RNLV,RTNLV,RNSO=gecros_utils.rnacc(FNSH,NUPT,RWST,STEMNC,LNCMIN,RNCMIN,LNC,RNC,NLV,     \
              NRT,WLV,WRT,DELT,CB,CX,TM,DS,SEEDNC,RWSO,LNLV,LNRT)
    
        #*** Estimation of total seed number, and 1000-seed weight
        RNRES  = NUPT-(LNCMIN*(RCLV+LCLV)+RNCMIN*(RCSRT+LCRT)+STEMNC*RCSST)/CFV
        RNREOE = gecros_utils.insw (DS-ESD, RNRES, 0.)
        RNREOF = gecros_utils.insw (DS-1.0, RNRES, 0.)
                
        #*** Plant height or stem length (m)
        RDCDTP = (DCDTC-DCDTP)/DELT
        
    
        RHT    = min(HTMX-HT, FDH*HTMX*IFSH)
        
        #*** Rooting depth (cm)
        KR     = -math.log(0.05)/RDMX 
        RRD    = gecros_utils.insw(RD-RDMX, min((RDMX-RD)/DELT,(RWRT+LWRT)/(WRB+KR* \
                    (WRT+WRTD))), 0.)
        
        #*** Soil water availability and water balance
        WCLL   = min(WCMAX, (WLL+WCMIN*10.*(RDMX1-RD))/10./(RDMX1-RD))
    
        RRUL   = min(10.*(WCMAX-WCUL)*      RD /TCP, RFIR)
        RRLL   = min(10.*(WCMAX-WCLL)*(RDMX1-RD)/TCP, RFIR-RRUL)
    
        RWUL   = RRUL+10.*(WCLL-WCMIN)*RRD-gecros_utils.insw(WSWI,0.,ATCAN+AESOIL)+.1
        RWLL   = RRLL-10.*(WCLL-WCMIN)*RRD
        RWUG   = max (0., RFIR-RRUL-RRLL)   
        ETA    = ATCAN+AESOIL
        
        #*** Soil mineral nitrogen
        LEAUL  = max(0.,(NSUPN-NUPTN)*TCP             -RD /RDMX1*RN)\
                  *min((RFIR-RRUL)/WCMAX/RD/10.,1.)
        LEALL  = max(0.,NNLL+(MINNLL-DENIUL)*TCP-(RDMX1-RD)/RDMX1*RN)\
                  *min(RWUG/WCMAX/(RDMX1-RD)/10.,1.)
    
        VOLA   = gecros_utils.insw (RAIN-1., 0.15, 0.) * SFERNA
        
        RSFNA  = FERNA - SFERNA/3.
    
        LAYNA  = RRD/(RDMX1-RD)*NALL
        LAYNN  = RRD/(RDMX1-RD)*NNLL

        RNAUL  =FERNA+MINAUL       +LAYNA-gecros_utils.insw(NSWI,0.,NUPTA)-NITRUL-VOLA
        RNALL  =      MINALL       -LAYNA                    -NITRLL
        RNNUL  =FERNN+MINNUL+NITRUL+LAYNN-gecros_utils.insw(NSWI,0.,NUPTN)-DENIUL-LEAUL
        RNNLL  =LEAUL+MINNLL+NITRLL-LAYNN                    -DENILL-LEALL
        
        #*** Maintenance and total respiration (g CO2 m-2 d-1)
        RMUN   = 44./12.*2.05*NUPTN
        RMUA   = 44./12.*0.17*NUPTA
        RMUS   = 0.06* 0.05/0.454*YGV*ASSA
        RMLD   = 0.06*(1.-FCSH)*ASSA
        RRMUL  = (RMUN+RMUA+RMUS+RMLD - RMUL)/DELT                   
                
        #*** Leaf area development
        SLA    = LAI /WLV   
        RLAI=gecros_utils.rlaic(DS,SLA0,RWLV,LAIC,KN,NLV,RNLV,SLNB,RSLNB)

        #*** Amount of seed protein
        PSO    = 6.25*WSO*ONC   
        PSOPC  = 6.25*ONC*100.
        #*** Daily carbon flow for seed filling
        RDCDSR = max(0., (DCDSC-RCSO/YGO))-(FLWCS-min(DCDSC,DCSS))
        
        #*** Daily carbon flow for structural stem growth
        RDCDTR = max(0., (DCDTC-RCSST/YGV))-(FLWCT-min(DCDTC,DCST))

        #*** Crop carbon balance check
        CCHKIN = CTOT + CLVD+CRTD -CLVI-CRTI
        CCHK   = (CCHKIN-((TPCAN-TRESP)*12./44.))/gecros_utils.notnul(CCHKIN)*100.
    
        #*** Crop nitrogen balance check
        NCHKIN = NTOT + NLVD+NRTD -NLVI-NRTI
        NCHK   = (NCHKIN-TNUPT)/gecros_utils.notnul(TNUPT)*100.
        
        #INTEGRATION #############################################################################
        #TRANSLATION_GENERAL DRIVER='EUDRIV' #simple rectangular (euler) integration
        TIME+=1
        #*** Carbon accumulation
        CLV    += RCLV#intgrl (CLVI, RCLV ,CLV)
        CLVD   += LCLV#intgrl (ZERO, LCLV ,CLVD)
        CSST   += RCSST#intgrl (ZERO, RCSST,CSST)
        CSO    += RCSO#intgrl (ZERO, RCSO ,CSO)
        CSRT   += RCSRT#intgrl (CRTI, RCSRT,CSRT)
        CRTD   += LCRT#intgrl (ZERO, LCRT ,CRTD)
        CLVDS  += LVDS#intgrl (ZERO, LVDS ,CLVDS)
        
        #*** Nitrogen accumulation
        NRT    += RNRT#intgrl (NRTI, RNRT ,NRT)
        NST    += RNST#intgrl (ZERO, RNST ,NST)
        NLV    += RNLV#intgrl (NLVI, RNLV ,NLV)
        NSO    += RNSO#intgrl (ZERO, RNSO ,NSO)
        TNLV   += RTNLV#intgrl (NLVI, RTNLV,TNLV)
        NLVD   += LNLV#intgrl (ZERO, LNLV ,NLVD)
        NRTD   += LNRT#intgrl (ZERO, LNRT ,NRTD)

        SLNB   += RSLNB#intgrl (SLNBI, RSLNB,SLNB)
        
        #*** Estimation of total seed number, and 1000-seed weight
        NREOE  += RNREOE#intgrl (ZERO, RNREOE,NREOE)
        NREOF  += RNREOF#intgrl (ZERO, RNREOF,NREOF)
        
        #*** Daily carbon flow for seed filling
        DCDSR  += RDCDSR#intgrl(ZERO, RDCDSR,DCDSR)

        #*** Daily carbon flow for structural stem growth
        DCDTR  += RDCDTR#intgrl(ZERO, RDCDTR,DCDTR)

    
        #*** Dynamics of carbon-reserve pool in stems and roots
        CRVS   += RCRVS#intgrl (ZERO, RCRVS,CRVS)
        CRVR   += RCRVR#intgrl (ZERO, RCRVR,CRVR)          

        #*** Leaf area development
        LAIC   += RLAI#intgrl(LAII, RLAI,LAIC)
        
        #*** Nitrogen fixation (g N m-2 d-1)
        NDEMP  += RNDEMP#intgrl (ZERO, RNDEMP,NDEMP)
        NSUPP  += RNSUPP#intgrl (ZERO, RNSUPP,NSUPP)
        NFIXT  += NFIX#intgrl (ZERO, NFIX,NFIXT)
        NFIXR  += RNFIXR#intgrl (ZERO, RNFIXR,NFIXR)
        
        RMUL   += RRMUL#intgrl(ZERO, RRMUL,RMUL) 
              
        DCDTP  += RDCDTP#intgrl(ZERO, RDCDTP,DCDTP)
        RD     += RRD#intgrl (RDI, RRD,RD)
        HT     += RHT#intgrl (HTI, RHT,HT) 
        
        #*** Canopy photosynthesis, transpiration and soil evaporation
        TDAPAR += DAPAR#intgrl (ZERO, DAPAR ,TDAPAR)
        TPCAN  += APCAN#intgrl (ZERO, APCAN ,TPCAN)
        TRESP  += RESTOT#intgrl (ZERO, RESTOT,TRESP)
        TTCAN  += ATCAN#intgrl (ZERO, ATCAN ,TTCAN)
        TNUPT  += NUPT#intgrl (ZERO, NUPT  ,TNUPT)
        
        #********************** THE EXAMPLE SOIL MODEL *************************
        #*** Soil temperature
        TSOIL  += RTSOIL#intgrl (TSOILI, RTSOIL,TSOIL)
        
        #*** Soil water availability and water balance
        WUL    += RWUL#intgrl (WULI, RWUL,WUL)
        WLL    += RWLL#intgrl (WLLI, RWLL,WLL)
        CETA   +=ETA
        WUG    +=RWUG
        
        #*** Soil organic carbon
        DPM    += RDPM#intgrl (DPMI, RDPM,DPM)
        RPM    += RRPM#intgrl (RPMI, RRPM,RPM)
        BIO    += RBIO#intgrl (BIOI, RBIO,BIO)
        HUM    += RHUM#intgrl (HUMI, RHUM,HUM)
    
        LITNT  += LITN#intgrl(ZERO, LITN,LITNT)
        
        #*** Soil organic nitrogen
        DPN    += RDPN#intgrl (DPNI, RDPN,DPN)
        RPN    += RRPN#intgrl (RPNI, RRPN,RPN)           
        
        SFERNA += RSFNA#intgrl (ZERO, RSFNA,SFERNA)

        TNLEA  += LEALL#intgrl (ZERO, LEALL,TNLEA)
        NAUL   += RNAUL#intgrl (NAULI, RNAUL,NAUL)
        NALL   += RNALL#intgrl (NALLI, RNALL,NALL)
        NNUL   += RNNUL#intgrl (NNULI, RNNUL,NNUL)
        NNLL   += RNNLL#intgrl (NNLLI, RNNLL,NNLL)
        TMDN   += MDNUL
        
        #*** Developmental stage (DS) & cumulative thermal units (CTDU)
        DS     += DVR#intgrl (ZERO, DVR,DS)
        CTDU   += TDU#intgrl (ZERO, TDU,CTDU)
        CVER   += VDU
        
        #OUTPUT ##################################################################################
        #print ('NTOT Nitrogen in living shoot and root ',NTOT,'gNm-2\n NLVD Nitrogen in dead leaves',NLVD,'gNm-2\n NRTD Nitrogen in dead roots',NRTD,'gNm-2\n NLVI Initial value of NLVI',NLVI,'gNm-2\n NRTI initial value of NRT',NRTI,'gNm-2 ')
        #print (('NUPTA Ammonium-N uptake by crop {0} gNm-2d-1  \n NUPTN Nitrate-N uptake by crop {1} gNm-2d-1 \n NDEM crop nitrogen demmand {2} gNm-2d-1'.format(NUPTA,NUPTN, NDEM)))
        #print (('NSUPA {0} NSUPN {1}'.format(NSUPA,NSUPN)))
        
        #print (('DFS Days from start {0}\n,DS {1}\n,carbon check CCHK {2}%\n,NCHK {3}%\n,HI harvest index {4}gg-1\n,WSO dry weight storage {5} gm-2\n,WSH dry w of shoot {6}gm-2\n,PSO prot content storage org {7}\n,TNUPT tot N uptake during growth {8} gNm-2\n,APCAN actual gross canopy photosynthesis {9} gCO2 m-2 d-1\n,PPCAN Potential canopy CO 2 assimilation {10}gCO2 m-2d-1\n,ATCAN Actual canopy transpiration {11}mm d-1\n,NUPT Crop N uptake {12}g N m-2 d-1\n,TSN total seed number {13} seeds m-2\n,ONC {14}\n,FCSH fraction of new carbon to shoot {15}gC gC-1\n,FNSH fraction of new absorved N partitioned to shoot {16}gN gN-1\n,LAI {17}'.format(DFS,DS,CCHK,NCHK,HI,WSO,WSH,PSO,TNUPT,APCAN,PPCAN,ATCAN,NUPT,TSN,ONC,FCSH,FNSH,LAI)))
        #Accumulate daily values in recorders         
        #k=0 # so all lists have equal lenght 
        #keys=dir(gecros)
        #for i in range(len(keys)):
            #if '_' in keys[i] or 'run' in keys[i] or 'd' in keys[i] or 'bristo' in keys[i]:# or keys[i] in ['i','run','bristo','values','keys','dict2','self','run','weather','afgen','VLS','IRRIG','FERTNA','FERTNN']:
                #continue
            #k=float(repr(eval('gecros.{0}'.format(keys[i]))))
            #if keys[i] not in self.d.keys():
                #self.d[keys[i]]=[k]
            #else:self.d[keys[i]].append(k)
        dict2=locals()
        for i in dict2.keys():
            if i in ['WEATHER','d','i','run','bristo','pheno','vern','values','keys','dict2','self','run','weather','afgen','VLS','IRRIG','FERTNA','FERTNN','fenologia_soja','k','t']:continue 
            
            k=float(dict2[i])
            if i not in d.keys():
                d[i]=[k]
            else:d[i].append(k)        
    print ('fisiological maturity ',DFS, datetime.datetime.fromordinal(TIME).isoformat())#DateTime.DateTimeFromJDN(TIME))
    #END
    #write output files 
    #self.keys=numpy.array(self.d.keys())
    #print (self.d.values())
    #self.d=numpy.array(self.d.values(),'f').T
    #fid=open('output_modelo_trigo.txt','w')
    #self.keys.tofile(fid, sep=",") 
    #fid.write('\n')
    #numpy.savetxt(fid, self.d, fmt="%5.5f",delimiter=',')  
    #fid.close()
    return d
    
    
def test(**kwargs):
    print (kwargs)
    d=run(**kwargs)
    pylab.plot(d['TIME'],d['DS'],d['TIME'],d['LAI'],d['TIME'],d['DLAI']);pylab.xlabel('DOY');pylab.ylabel('LAI,DLAI, DS')
    pylab.figure()
    #pylab.plot(d['TIME'],numpy.array(d['WSO'])*10.);pylab.xlabel('DOY');pylab.ylabel('Grain yield (Kg Ha -1)')
    #pylab.figure()
    pylab.plot(d['TIME'],d['TNUPT'], d['TIME'],(numpy.array(d['NNUL'])+numpy.array(d['NNLL']))/(d['RDMX'][-1]*.01*1.4), d['TIME'],(numpy.array(d['NAUL'])+numpy.array(d['NALL']))/(d['RDMX'][-1]*.01*1.4)) ;  pylab.xlabel('DOY');pylab.ylabel('TNUPT, NN_ppm, NA_ppm')
    pylab.figure()
    pylab.plot(d['TIME'],d['WSHH'],d['TIME'],d['WSH'],d['TIME'],d['WLV'],d['TIME'],d['WLVD'],d['TIME'],d['WSO']) 
    pylab.xlabel('DOY');pylab.ylabel('WSHH,WSH,WLV,WLVD,WSO')
    


if __name__ == '__main__' :
    if len(sys.argv)==2 : 
        scenarioName = sys.argv[1]
        print ('Running scenario %s' %(scenarioName))
        run()
    else :
        print ('Please re-enter your request.')
        print ('enter, for example')
        print ('python modelo_inia_trigo.py scenario_name')

 
        
#        PLANT
#        =====
#        CLVI  -Initial value, amount of carbon in the living leaves               gCm-2
#        CRTI  -Initial value, amount of carbon in living roots                    gCm-2
#        NRTI  -Initial value, nitrogen in living roots (Ng)                       gNm-2
#        NLVI  -Initial value, nitrogen in living leaves(Nlv)                      gNm-2
#        SLNBI -SLN in bottom leaves of the canopy (nbot)                          gNm-2 
#        LAII  -Initial value, of green LAI                                        m2m-2
#        HTI   -Initial value, of plant height                                     m
#        RDI   -Initial value, of rooting depth to the soil                        cm
#        
#        SOIL
#        ====
#        TSOILI-Initial value, of soil temperature                                 C
#        WULI  -Initial value, water content in the upper soil layer               mm
#        WLLI  -Initial value, water content in the lower soil layer               mm
#        DPMI  -Initial value, decomposable plant material                         gCm-2
#        RPMI  -Initial value, resistant plant material (difficult to decompose)   gCm-2
#        BIOI  -Initial value, microbial biomass in the soil                       gCm-2
#        HUMI  -Initial value, humidified organic matter in the soil               gCm-2
#        DPNI  -Initial value, organic N in decomposable plant material            gNm-2
#        RPNI  -Initial value, organic N in resistant plant material               gNm-2
#        NAULI -Initial value, amonium N in the upper soil layer                   gNm-2
#        NALLI -Initial value, amonium N in the lower soil layer                   gNm-2
#        NNULI -Initial value, nitrate N in the upper soil layer                   gNm-2
#        NNLLI -Initial vrunalue, nitrate N in the lower soil layer                   gNm-2       

#        PESOIL  #Potential soil evaporation                                       mmd-1
#        WSO     #Dry weight of storage organs                                     gm-2
#        WSHH    #Dry weight of shoot organs (excludding shedded leaves)           gm-2
#        WLVD    #Dry weight of dead leaves                                        gm-2
#        RCLV    Rate of change in CLV                                             gCm-2d-1
#        APCAN   #Actual gross canopy photosynthesis                               gCO2m-2d-1
#        ASSA    #Assimilates available from current photosynthesis for growth     gCO2m-2d-1
#        FCLV    #Fraction of new shoot carbon partitioned to leaves               gCg-1C
#        NCR     #Intermediate variable (nu)                                       gNg-1C
#        DERI    #First order derivative of SHSA with respect to crop N/C ratio    gCg-1Nd-1
#        SHSA    #Relative shoot activity                                          gCg-1Cd-1
#        FCSH    #Fraction of new carbon partitioned to shoot                      gCg-1C
#        DCSS    #Daily carbon supply from current photosynthesis for shoot growth gCm-2d-1
#        CSRTN   #Nitrogen determined CSRT                                         gCm-2
#        FNSH    #Fraction of newly absorbed N partitioned to the shoot            gNg-1N
#        KCRN    #Extinction coefficient of root nitrogen                          m2g-1C
#        LAI     #Green leaf area index                                            m2leafm-2
#        KN      #Leaf nitrogen extinction coefficient in the canopy               m2m-2leaf
#        TLAI    #Total leaf area index (green and seneced)                        m2leafm-2
#        #RLAI   #Rate of change in LAIC                                           m2leafm-2d-1
#        NSUP    #Nitrogen supply to crop                                          gNm-2d-1
#        DWSUP   #Daily water supply for evapotranspiration                        mm d-1
#        WCUL    #Soil water content in the upper soil layer                       m3m-3

