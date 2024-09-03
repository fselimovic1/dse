from constZ2nd import constZ2nd
from constZ2ndS import constZ2ndS

def td(settings, ppc):
    
    if settings["tdmethod"] == "constZ2nd": 
        results = constZ2nd(settings, ppc);
    elif settings["tdmethod"] == "constZ2ndS":
        results = constZ2ndS(settings, ppc);
    return results;
