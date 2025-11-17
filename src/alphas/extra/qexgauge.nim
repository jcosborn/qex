import gauge/[gaugeAction, gaugeUtils]

import alphasnorm

proc gaugeActionOneLoopHISQ*(c: GaugeActionCoeffs; g: array|seq): auto =
  tic("gaugeActionOneLoopHISQ")
  const nc = g[0][0].nrows
  let lo = g[0].l
  let nd = lo.nDim
  let t = newTransporters(g, g[0], 1)
  let t2 = newTransporters(g, g[0], 1)
  let td = newTransporters(g, g[0], -1)
  var pl = 0.0
  var rt = 0.0
  var pg = 0.0
  toc("gaugeActionOneLoopHISQ setup")
  threads:
    tic()
    toc("gaugeActionOneLoopHISQ zero")
    for mu in 1..<nd:
      for nu in 0..<mu:
        tic()
        var tpl = redot2(t[mu]^*g[nu], t[nu]^*g[mu])
        if threadNum == 0: pl += tpl
        toc("gaugeActionOneLoopHISQ pl")
        if c.rect != 0:
          var tr1 = redot2(t[mu]^*t[nu]^*g[nu], t2[nu]^*t[nu]^*g[mu])
          var tr2 = redot2(t2[mu]^*t[mu]^*g[nu], t[nu]^*t[mu]^*g[mu])
          if threadNum == 0: rt += tr1 + tr2
          toc("gaugeActionOneLoopHISQ rt")
        if c.pgm != 0:
          for sg in 0..<nu:
            var ts1 = redot2(t[mu]^*t[nu]^*g[sg], t[sg]^*t[nu]^*g[mu])
            var ts2 = redot2(t[mu]^*t[sg]^*g[nu], t[nu]^*t[sg]^*g[mu])
            var ts3 = redot2(t[nu]^*t[mu]^*g[sg], t[sg]^*t[mu]^*g[nu])
            var ts4 = redot2(t[nu]^*t[sg]^*g[mu], t[mu]^*t[sg]^*g[nu])
            var ts5 = redot2(t[sg]^*t[mu]^*g[nu], t[nu]^*t[mu]^*g[sg])
            var ts6 = redot2(t[sg]^*t[nu]^*g[mu], t[mu]^*t[nu]^*g[sg])
            var ts7 = redot2(t[mu]^*td[nu]^*g[sg], t[sg]^*td[nu]^*g[mu])
            var ts8 = redot2(t[mu]^*td[sg]^*g[nu], t[nu]^*td[sg]^*g[mu])
            if threadNum == 0:
              pg += ts1 + ts2 + ts3 + ts4 + ts5 + ts6 + ts7 + ts8
          toc("gaugeActionOneLoopHISQ pg")
    toc("gaugeActionOneLoopHISQ work")
  toc("gaugeActionOneLoopHISQ threads")
  result = (1.0/nc.float) * (c.plaq*pl + c.rect*rt + c.pgm*pg)