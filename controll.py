import biglimpcopy1.petprocessing as pp
import numpy as np
import json
import requests

__author__="noahtingb"

def getdata(cor=[57,12],types=["temperature_2m"],forcastdays=14):
        return requests.get('https://api.open-meteo.com/v1/forecast?latitude='+str(cor[0])+'&longitude='+str(cor[1])+'&hourly='+",".join(types)+"&forecast_days="+str(forcastdays))
def gethourlyparams1():        
    return ["temperature_2m", "relative_humidity_2m", "cloud_cover", "cloud_cover_low", "cloud_cover_mid", "cloud_cover_high", "wind_speed_10m", "wind_gusts_10m", "shortwave_radiation", "direct_radiation", "diffuse_radiation"]

def dumpa(name,result):
    with open(name, 'w') as f:
            json.dump(result,f)
def loada(name):
    with open(name,"r") as f:
        return json.load(f)

def checkwater(lat,lon):
    return requests.get("https://is-on-water.balbona.me/api/v1/get/"+str(lat)+"/"+str(lon)).json()
def createrandomcordinates(n=200):
    lats,lons=np.random.uniform(-90+1e-10,90,n),np.random.uniform(-180+1e-10,180,n)
    watercords,landcords=[],[]
    for i in range(n):
        print(i/n,i)
        response=checkwater(lats[i],lons[i])
        print(response)
        if response["isWater"]:
            watercords.append([lats[i],lons[i]])
        else:
            landcords.append([lats[i],lons[i]])
    dumpa("landandwatercord1000.json",{"landcords":landcords,"watercords":watercords})
def combine():
    listan=loada("landandwatercord1000.json")["landcords"]+loada("landandwatercord.json")["landcords"]
    dumpa("combineisland.json",{"cord":listan})
def createmuchdata(vikt=25):
    cordaa=loada("combineisland.json")["cord"]
    cords=cordaa[vikt:min(vikt+3,320)]
    for i in range(len(cords)):
        data=getdata(cords[i],gethourlyparams1()).json()
        dumpa(f"weatherdata//weatherdata{vikt+i}.json",{"weatherdata":data,"cords":cords[i]})
        print((vikt+i)/len(cordaa),i+vikt,len(cordaa))
    print("made")
    #dumpa("weatherdata.json",{"weatherdata":data,"cords":cords})
def createparams(data,time=[2025,4,2,0],cords=[57,12]):
    time[3]=np.random.randint(0,23)
    time[2]=np.random.randint(2,15)
    timestr=list(str(time[0])+"-04-02T00:00")
    for i in range(len(str(time[1]))):
        timestr[6-i]=str(time[1])[len(str(time[1]))-1-i]
    for i in range(len(str(time[2]))):
        timestr[9-i]=str(time[2])[len(str(time[2]))-1-i]
    for i in range(len(str(time[3]))):
        timestr[12-i]=str(time[3])[len(str(time[3]))-1-i]
    timestr="".join(timestr)
    if timestr in data["hourly"]["time"]:
        indexforsearch=data["hourly"]["time"].index(timestr)
        params={"year":time[0],"month":time[1],"day":time[2],"hour":time[3],"Ta":data["hourly"]["temperature_2m"][indexforsearch],"RH":data["hourly"]["relative_humidity_2m"][indexforsearch],"Ws":data["hourly"]["wind_speed_10m"][indexforsearch],"sky":"Clear (100%)","c1":data["hourly"]["cloud_cover_low"][indexforsearch],"c2":data["hourly"]["cloud_cover_mid"][indexforsearch],"c3":data["hourly"]["cloud_cover_high"][indexforsearch], "loc":{'latitude': cords[0],'longitude': cords[1], 'altitude': data["elevation"]},"UTC":0}
        return params
    return None
def controll(time=[2025,4,2,20]):
    data=loada(f"weatherdata//weatherdata{322}.json")
    cords=data["cords"]
    data=data["weatherdata"]
    params=createparams(data,time,cords)
    #print(params)
    petresult=pp.indexflask(params)
    return petresult,params,data