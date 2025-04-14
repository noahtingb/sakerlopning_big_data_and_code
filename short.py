def Lvikt_veg():
        viktonlywall=(4.4897-(63.227*0.6**6-161.51*0.6**5+156.91*0.6**4-70.424*0.6**3+16.773*0.6**2-0.4863*0.6))/4.4897
        viktaveg=(4.4897-(63.227*1**6-161.51*1**5+156.91*1**4-70.424*1**3+16.773*1**2-0.4863*1))/4.4897
        viktwall=viktonlywall-viktaveg
        svfvegbu=(1+0.6-1)  # Vegetation plus buildings
        viktsky=(63.227*svfvegbu**6-161.51*svfvegbu**5+156.91*svfvegbu**4-70.424*svfvegbu**3+16.773*svfvegbu**2-0.4863*svfvegbu)/4.4897
        viktrefl=(4.4897-(63.227*svfvegbu**6-161.51*svfvegbu**5+156.91*svfvegbu**4-70.424*svfvegbu**3+16.773*svfvegbu**2-0.4863*svfvegbu))/4.4897
        viktveg=(4.4897-(63.227*svfvegbu**6-161.51*svfvegbu**5+156.91*svfvegbu**4-70.424*svfvegbu**3+16.773*svfvegbu**2-0.4863*svfvegbu))/4.4897
        viktveg=viktveg-viktwall
        return viktveg*1e10,viktwall*1e10,viktsky*1e10,viktrefl*1e10
def Kvikt_veg():
    # Least
    viktwall=(4.4897-(63.227*0.6**6-161.51*0.6**5+156.91*0.6**4-70.424*0.6**3+16.773*0.6**2-0.4863*0.6))/4.4897
    
    svfvegbu=(1+0.6-1)  # Vegetation plus buildings
    viktveg=(4.4897-(63.227*svfvegbu**6-161.51*svfvegbu**5+156.91*svfvegbu**4-70.424*svfvegbu**3+16.773*svfvegbu**2-0.4863*svfvegbu))/4.4897
    viktveg=viktveg-viktwall    
    return viktveg,viktwall
def Kvikt_veg1(svf=0.6,svfveg=1,vikttot=4.4897):
    # Least
    viktwall=(vikttot-(63.227*svf**6-161.51*svf**5+156.91*svf**4-70.424*svf**3+16.773*svf**2-0.4863*svf))/vikttot
    
    svfvegbu=(svfveg+svf-1)  # Vegetation plus buildings
    viktveg=(vikttot-(63.227*svfvegbu**6-161.51*svfvegbu**5+156.91*svfvegbu**4-70.424*svfvegbu**3+16.773*svfvegbu**2-0.4863*svfvegbu))/vikttot
    viktveg=viktveg-viktwall    
    return viktveg,viktwall
print(Kvikt_veg(),Kvikt_veg1(),Lvikt_veg())