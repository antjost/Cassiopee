# - enforceh (array) -
import Geom as D
import Converter as C

a = D.line((0,0,0), (1,0,0), N=30)
D.setH(a, 0, 0.01); D.setH(a, -1, 0.1)

b = D.enforceh(a, N=40)
C.convertArrays2File(b, 'out.plt')
