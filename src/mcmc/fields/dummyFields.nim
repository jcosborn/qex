import ../mcmcTypes

proc newDummyField(l: Layout; info: JsonNode): auto = 
    result = l.newLatticeField(info)
    result.du = l.newGauge()
    result.su = l.ColorVector()