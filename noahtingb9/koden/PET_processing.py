'''A lower controll function for the PET Calculation'''
import noahtingb9.koden.Solweig as so
import noahtingb9.koden.PET_calculations as p
__author__="noahtingb"

def indexflask(form):
        month = int(form["month"])
        day = int(form["day"])
        hour = int(form["hour"])
        year = int(form["year"])
        minut = 30
        Ta = float(form["Ta"])        
        RH = float(form["RH"])
        Ws = float(form["Ws"])
        location = form["loc"]

        if month > 12 or month < 0:
            print("petresult.html","Incorrect month filled in")
        if day > 31 or day < 0:
            print("petresult.html","Incorrect day filled in")
        if hour > 23 or hour < 0:
            print("petresult.html","Incorrect hour filled in")
        if Ta > 60 or Ta < -75:
            print("petresult.html", "Unreasonable air temperature filled in",Ta)
        if RH > 100 or RH < 0:
            print("petresult.html", "Unreasonable relative humidity filled in")
        if Ws > 100 or Ws < 0:
            print("petresult.html", "Unreasonable Wind speed filled in")
 
        # Main calculation
        if Ta is not None and RH is not None and Ws is not None:
            Tmrt, resultPET = petcalc(Ta, RH, Ws, year, month, day, hour,location)
        return resultPET,Tmrt

def petcalc(Ta, RH, Ws, year, month, day, hour,location):    

    Fside,Fup,Fcyl = 0.22,0.06,0.28 #Ståendes. Vid Liggande:    Fside,Fup,Fcyl = 0.166666, 0.166666, 0.2
    
    Tmrt = so.Solweig1D_2020a_calc(Fside, Fup, Fcyl,location,Ta, RH, year, month, day, hour,minu=30)

    WsPET = (1.1 / 10) ** 0.2 * Ws #corretion from 10 meters height to 1.1 meters height 
    
    mbody, ht, clo, age, activity, sex = 75., 1.8, 0.9, 35, 80.,  1#[kg], [m], [1], [years], [W], [m 1/f 2]
    
    resultPET = p._PET(Ta, RH, Tmrt, WsPET, mbody, age, ht, activity, clo, sex) #get Pet

    return Tmrt, resultPET
