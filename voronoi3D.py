import pyvoro
from pyvoro import voroplusplus
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from mpl_toolkits.mplot3d import Axes3D
import scipy as sp
import scipy.spatial as sptl
import scipy.sparse as sprs
import numpy as np
import pickle

### COORD 3D RECTANGLE
x=np.arange(-24, 24, 3)
y=np.arange(-24, 24, 3)
#print(y)
xx,yy=np.meshgrid(x,y)
#print(yy)
xx=(np.sqrt(3)/ 2)*xx
#print('xx',xx)
a1=np.array([0,1])
yy=yy+np.matlib.repmat(a1, np.size(xx,0), int(np.size(xx,0)/2))

A= np.asarray(xx.T).reshape(-1)
#print('A',A)
B=np.asarray(yy.T).reshape(-1)
#for i in range(len(C)):
#    if i%2==0:
#        C[i]=C[i]*0.8
#print(B)
C=np.zeros(len(B))
C1=np.ones(len(B))*3
C2=C1*2
C3=C1*3
C4=C1*4
C5=C1*5
C6=C1*6
M=np.column_stack((A,B))
U=np.column_stack((M,C1))
U=U.tolist()
J=np.column_stack((M,C2))
J=J.tolist()
W=np.column_stack((M,C3))
W=W.tolist()
Z=np.column_stack((M,C4))
Z=Z.tolist()
S=np.column_stack((M,C5))
S=S.tolist()
L=np.column_stack((M,C6))
L=L.tolist()
Q=np.column_stack((M,C))
Q=Q.tolist()
for i in U:
    Q.append(i)
for i in J:
    Q.append(i)
for i in W:
    Q.append(i)
for i in Z:
    Q.append(i)
for i in S:
    Q.append(i)
for i in L:
    Q.append(i)
##################
shell_cone=[]
volume_cone=[]
x_cone=np.arange(-5,5,0.5)
y_cone=np.arange(-5,5,0.5)
z_cone=np.arange(-10,0,0.5)
for i in x_cone:
    for j in y_cone:
        for k in z_cone:
            if -0.1<i**2+j**2-k**2<0.1:
                shell_cone.append([i,j,k])
            if -5<i**2+j**2-k**2<5:
                volume_cone.append([i,j,k])
#coord_cone=[[-1.2184253819214064e-06, 1.947204351425171, 2.595390796661377], [-1.2184253819214064e-06, 3.8195760250091553, 2.0274128913879395], [-1.2184253819214064e-06, 5.545164585113525, 1.1050671339035034], [-1.2184253819214064e-06, 7.05765438079834, -0.13620156049728394], [-1.2184253819214064e-06, 8.298924446105957, -1.6486923694610596], [0.3798796534538269, 1.9097893238067627, 2.595390796661377], [0.7451613545417786, 3.7461841106414795, 2.0274128913879395], [1.0818067789077759, 5.438615322113037, 1.1050671339035034], [1.3768787384033203, 6.922043323516846, -0.13620156049728394], [1.6190383434295654, 8.1394624710083, -1.6486923694610596], [0.7451618909835815, 1.7989823818206787, 2.595390796661377], [1.4616873264312744, 3.5288283824920654, 2.0274128913879395], [2.1220414638519287, 5.123063564300537, 1.1050671339035034], [2.7008464336395264, 6.520421981811523, -0.13620156049728394], [3.175859212875366, 7.667204856872559, -1.6486923694610596], [1.0818079710006714, 1.6190412044525146, 2.595390796661377], [2.122041940689087, 3.1758615970611572, 2.0274128913879395], [3.0807268619537354, 4.610634803771973, 1.1050671339035034], [3.9210214614868164, 5.868224620819092, -0.13620156049728394], [4.610633373260498, 6.900301933288574, -1.6486923694610596], [1.3768806457519531, 1.3768813610076904, 2.595390796661377], [2.7008471488952637, 2.70084810256958, 2.0274128913879395], [3.9210219383239746, 3.921022415161133, 1.1050671339035034], [4.990513801574707, 4.990515232086182, -0.13620156049728394], [5.868223667144775, 5.868224143981934, -1.6486923694610596], [1.6190403699874878, 1.0818085670471191, 2.595390796661377], [3.175860643386841, 2.1220428943634033, 2.0274128913879395], [4.610633850097656, 3.0807273387908936, 1.1050671339035034], [5.868223667144775, 3.9210216999053955, -0.13620156049728394], [6.900300979614258, 4.610633850097656, -1.6486923694610596], [1.7989813089370728, 0.7451623678207397, 2.595390796661377], [3.528827428817749, 1.4616882801055908, 2.0274128913879395], [5.123062610626221, 2.122042179107666, 1.1050671339035034], [6.520421504974365, 2.7008461952209473, -0.13620156049728394], [7.667203426361084, 3.175859212875366, -1.6486923694610596], [1.9097883701324463, 0.3798801004886627, 2.595390796661377], [3.746182918548584, 0.7451618909835815, 2.0274128913879395], [5.4386138916015625, 1.0818073749542236, 1.1050671339035034], [6.922041893005371, 1.3768794536590576, -0.13620156049728394], [8.139459609985352, 1.619038701057434, -1.6486923694610596], [1.9472031593322754, -6.48571358397021e-07, 2.595390796661377], [3.8195748329162598, -5.742069788539084e-07, 2.0274128913879395], [5.545162677764893, -4.998424856239581e-07, 1.1050671339035034], [7.057652473449707, -5.742069788539084e-07, -0.13620156049728394], [8.298920631408691, -7.229358516269713e-07, -1.6486923694610596], [1.9097881317138672, -0.37988144159317017, 2.595390796661377], [3.746182680130005, -0.7451629042625427, 2.0274128913879395], [5.4386138916015625, -1.0818082094192505, 1.1050671339035034], [6.9220404624938965, -1.376880168914795, -0.13620156049728394], [8.139458656311035, -1.6190398931503296, -1.6486923694610596], [1.798980951309204, -0.7451635003089905, 2.595390796661377], [3.5288267135620117, -1.4616888761520386, 2.0274128913879395], [5.1230621337890625, -2.122042655944824, 1.1050671339035034], [6.520420074462891, -2.7008473873138428, -0.13620156049728394], [7.667201995849609, -3.175860643386841, -1.6486923694610596], [1.61903977394104, -1.081809401512146, 2.595390796661377], [3.1758596897125244, -2.1220428943634033, 2.0274128913879395], [4.610633373260498, -3.08072829246521, 1.1050671339035034], [5.868223190307617, -3.921022415161133, -0.13620156049728394], [6.900299549102783, -4.610634803771973, -1.6486923694610596], [1.3768796920776367, -1.3768819570541382, 2.595390796661377], [2.7008466720581055, -2.7008485794067383, 2.0274128913879395], [3.921020984649658, -3.921022891998291, 1.1050671339035034], [4.990512847900391, -4.990513801574707, -0.13620156049728394], [5.868221759796143, -5.868223667144775, -1.6486923694610596], [1.081807255744934, -1.6190414428710938, 2.595390796661377], [2.1220414638519287, -3.175861358642578, 2.0274128913879395], [3.080725908279419, -4.610634803771973, 1.1050671339035034], [3.921020030975342, -5.868223667144775, -0.13620156049728394], [4.610630989074707, -6.9003005027771, -1.6486923694610596], [0.7451614141464233, -1.7989823818206787, 2.595390796661377], [1.4616870880126953, -3.5288279056549072, 2.0274128913879395], [2.1220407485961914, -5.123062610626221, 1.1050671339035034], [2.7008447647094727, -6.520420074462891, -0.13620156049728394], [3.1758570671081543, -7.667202472686768, -1.6486923694610596], [0.3798793852329254, -1.9097893238067627, 2.595390796661377], [0.7451612949371338, -3.746183395385742, 2.0274128913879395], [1.0818066596984863, -5.4386138916015625, 1.1050671339035034], [1.3768779039382935, -6.922040939331055, -0.13620156049728394], [1.6190366744995117, -8.139457702636719, -1.6486923694610596], [-1.2904087043352774e-06, -1.9472042322158813, 2.595390796661377], [-8.070397825576947e-07, -3.819575309753418, 2.0274128913879395], [-9.557686553307576e-07, -5.545162677764893, 1.1050671339035034], [-1.476319994253572e-06, -7.057651042938232, -0.13620156049728394], [-2.0712354853458237e-06, -8.298918724060059, -1.6486923694610596], [-1.263082504272461, 6.807505130767822, -1.6797120571136475], [-1.0598063468933105, 5.785562038421631, -0.41007059812545776], [-0.8121118545532227, 4.540316104888916, 0.6318969130516052], [-0.5295193791389465, 3.1196258068084717, 1.4061485528945923], [-0.22288890182971954, 1.5780889987945557, 1.8829305171966553], [-2.5699353218078613, 6.411076545715332, -1.6797120571136475], [-0.3798818588256836, -1.9097890853881836, 2.595390796661377], [-0.745162844657898, -3.746182918548584, 2.0274128913879395], [-1.0818084478378296, -5.4386138916015625, 1.1050671339035034], [-1.3768812417984009, -6.9220404624938965, -0.13620156049728394], [-1.6190409660339355, -8.139456748962402, -1.6486923694610596], [-2.171193838119507, 5.448426723480225, -0.41007059812545776], [-1.6853232383728027, 4.275430679321289, 0.6318969130516052], [-1.130998134613037, 2.937169313430786, 1.4061485528945923], [-0.5295203328132629, 1.485073447227478, 1.8829305171966553], [-3.7743382453918457, 5.767310619354248, -1.6797120571136475], [-3.1954543590545654, 4.9009480476379395, -0.41007059812545776], [-2.490079641342163, 3.8452792167663574, 0.6318969130516052], [-1.6853234767913818, 2.640876293182373, 1.4061485528945923], [-0.8121134042739868, 1.334024429321289, 1.8829305171966553], [-4.830005645751953, 4.900947093963623, -1.6797120571136475], [-0.7451637983322144, -1.7989819049835205, 2.595390796661377], [-1.461688756942749, -3.528826951980591, 2.0274128913879395], [-2.122042655944824, -5.123062610626221, 1.1050671339035034], [-2.700847625732422, -6.520419597625732, -0.13620156049728394], [-3.175860643386841, -7.667199611663818, -1.6486923694610596], [-4.093225955963135, 4.16416597366333, -0.41007059812545776], [-3.1954545974731445, 3.2663936614990234, 0.6318969130516052], [-2.171194076538086, 2.2421324253082275, 1.4061485528945923], [-1.0598082542419434, 1.1307464838027954, 1.8829305171966553], [-5.6963701248168945, 3.8452792167663574, -1.6797120571136475], [-4.830008506774902, 3.2663943767547607, -0.41007059812545776], [-3.774341106414795, 2.561018943786621, 0.6318969130516052], [-2.5699381828308105, 1.7562625408172607, 1.4061485528945923], [-1.2630863189697266, 0.8830515146255493, 1.8829305171966553], [-6.340137004852295, 2.6408770084381104, -1.6797120571136475], [-1.0818095207214355, -1.6190409660339355, 2.595390796661377], [-2.122042655944824, -3.1758601665496826, 2.0274128913879395], [-3.080728054046631, -4.610634803771973, 1.1050671339035034], [-3.921022415161133, -5.868222713470459, -0.13620156049728394], [-4.610633850097656, -6.900297164916992, -1.6486923694610596], [-5.3774871826171875, 2.2421340942382812, -0.41007059812545776], [-4.204492092132568, 1.7562627792358398, 0.6318969130516052], [-2.866231918334961, 1.2019367218017578, 1.4061485528945923], [-1.414135456085205, 0.6004586219787598, 1.8829305171966553], [-6.736567497253418, 1.334024429321289, -1.6797120571136475], [-5.71462345123291, 1.1307467222213745, -0.41007059812545776], [-4.4693779945373535, 0.8830510973930359, 0.6318969130516052], [-3.0486881732940674, 0.6004583239555359, 1.4061485528945923], [-1.5071511268615723, 0.29382726550102234, 1.8829305171966553], [-6.870425701141357, -0.02505718357861042, -1.6797120571136475], [-3.5180989925720496e-06, -2.9403145163087174e-06, 2.7871735095977783], [-1.3768820762634277, -1.3768812417984009, 2.595390796661377], [-2.700847625732422, -2.7008469104766846, 2.0274128913879395], [-3.921022891998291, -3.921022653579712, 1.1050671339035034], [-4.990514278411865, -4.990512847900391, -0.13620156049728394], [-5.868222713470459, -5.86821985244751, -1.6486923694610596], [-5.828460693359375, -0.025057559832930565, -0.41007059812545776], [-4.558819770812988, -0.02505836822092533, 0.6318969130516052], [-3.1102964878082275, -0.02505836822092533, 1.4061485528945923], [-1.5385587215423584, -0.025058617815375328, 1.8829305171966553], [-6.736568450927734, -1.384138822555542, -1.6797120571136475], [-5.714623928070068, -1.1808619499206543, -0.41007059812545776], [-4.46937894821167, -0.9331681132316589, 0.6318969130516052], [-3.0486881732940674, -0.6505751013755798, 1.4061485528945923], [-1.5071513652801514, -0.34394460916519165, 1.8829305171966553], [-6.3401384353637695, -2.6909918785095215, -1.6797120571136475], [-1.6190416812896729, -1.0818085670471191, 2.595390796661377], [-3.175860643386841, -2.1220414638519287, 2.0274128913879395], [-4.610634803771973, -3.0807275772094727, 1.1050671339035034], [-5.868223667144775, -3.9210212230682373, -0.13620156049728394], [-6.900299072265625, -4.610630989074707, -1.6486923694610596], [-5.377488136291504, -2.2922494411468506, -0.41007059812545776], [-4.204493045806885, -1.8063796758651733, 0.6318969130516052], [-2.86623215675354, -1.2520536184310913, 1.4061485528945923], [-1.4141356945037842, -0.6505759358406067, 1.8829305171966553], [-5.6963725090026855, -3.895394802093506, -1.6797120571136475], [-4.830009460449219, -3.3165104389190674, -0.41007059812545776], [-3.7743422985076904, -2.6111364364624023, 0.6318969130516052], [-2.569938898086548, -1.8063794374465942, 1.4061485528945923], [-1.2630865573883057, -0.9331690073013306, 1.8829305171966553], [-4.830008506774902, -4.951062202453613, -1.6797120571136475], [-1.7989825010299683, -0.7451629042625427, 2.595390796661377], [-3.528826951980591, -1.461687445640564, 2.0274128913879395], [-5.1230621337890625, -2.122041702270508, 1.1050671339035034], [-6.520420551300049, -2.700845956802368, -0.13620156049728394], [-7.667200565338135, -3.1758573055267334, -1.6486923694610596], [-4.093227386474609, -4.214282035827637, -0.41007059812545776], [-3.195456027984619, -3.316511631011963, 0.6318969130516052], [-2.1711952686309814, -2.292250394821167, 1.4061485528945923], [-1.0598087310791016, -1.1808640956878662, 1.8829305171966553], [0.09599475562572479, -0.025060666725039482, 2.0439200401306152], [-3.774341344833374, -5.817426681518555, -1.6797120571136475], [-3.19545578956604, -4.951064586639404, -0.41007059812545776], [-2.4900810718536377, -3.895397901535034, 0.6318969130516052], [-1.6853244304656982, -2.6909942626953125, 1.4061485528945923], [-0.8121138215065002, -1.3841420412063599, 1.8829305171966553], [-1.9097895622253418, -0.37988102436065674, 2.595390796661377], [-3.7461819648742676, -0.7451618909835815, 2.0274128913879395], [-5.4386138916015625, -1.0818074941635132, 1.1050671339035034], [-6.922041416168213, -1.3768787384033203, -0.13620156049728394], [-8.139456748962402, -1.61903715133667, -1.6486923694610596], [-2.569938898086548, -6.461193561553955, -1.6797120571136475], [-2.1711952686309814, -5.4985432624816895, -0.41007059812545776], [-1.6853244304656982, -4.325549602508545, 0.6318969130516052], [-1.1309986114501953, -2.987287759780884, 1.4061485528945923], [-0.5295207500457764, -1.5351911783218384, 1.8829305171966553], [-1.2630860805511475, -6.857624053955078, -1.6797120571136475], [-1.0598080158233643, -5.835679054260254, -0.41007059812545776], [-0.8121129274368286, -4.590435028076172, 0.6318969130516052], [-0.5295199751853943, -3.1697447299957275, 1.4061485528945923], [-0.22288937866687775, -1.6282069683074951, 1.8829305171966553], [-1.947204351425171, -4.998424856239581e-07, 2.595390796661377], [-3.8195743560791016, -2.0238474007783225e-07, 2.0274128913879395], [-5.545162200927734, -2.0238474007783225e-07, 1.1050671339035034], [-7.057652473449707, 7.643529897904955e-07, -0.13620156049728394], [-8.298917770385742, 1.210539721796522e-06, -1.6486923694610596], [0.0959959626197815, -6.991482257843018, -1.6797120571136475], [0.09599646925926208, -5.949515342712402, -0.41007059812545776], [0.0959969013929367, -4.67987585067749, 0.6318969130516052], [0.09599703550338745, -3.2313530445098877, 1.4061485528945923], [0.09599661827087402, -1.6596145629882812, 1.8829305171966553], [1.455077886581421, -6.857624530792236, -1.6797120571136475], [1.251800775527954, -5.835679531097412, -0.41007059812545776], [1.0041069984436035, -4.590435028076172, 0.6318969130516052], [0.7215140461921692, -3.1697449684143066, 1.4061485528945923], [0.41488274931907654, -1.6282070875167847, 1.8829305171966553], [-1.9097893238067627, 0.3798799216747284, 2.595390796661377], [-3.7461819648742676, 0.7451614737510681, 2.0274128913879395], [-5.438612937927246, 1.081807017326355, 1.1050671339035034], [-6.922040939331055, 1.3768802881240845, -0.13620156049728394], [-8.139455795288086, 1.6190396547317505, -1.6486923694610596], [2.7619314193725586, -6.461195945739746, -1.6797120571136475], [2.3631882667541504, -5.498543739318848, -0.41007059812545776], [1.8773183822631836, -4.325549602508545, 0.6318969130516052], [1.3229928016662598, -2.9872887134552, 1.4061485528945923], [0.7215141654014587, -1.5351916551589966, 1.8829305171966553], [3.966334581375122, -5.817429542541504, -1.6797120571136475], [3.387449264526367, -4.951065540313721, -0.41007059812545776], [2.682074785232544, -3.895397901535034, 0.6318969130516052], [1.877319097518921, -2.690995216369629, 1.4061485528945923], [1.0041072368621826, -1.3841423988342285, 1.8829305171966553], [-1.7989823818206787, 0.7451618313789368, 2.595390796661377], [-3.5288267135620117, 1.4616870880126953, 2.0274128913879395], [-5.123061180114746, 2.1220409870147705, 1.1050671339035034], [-6.520419597625732, 2.7008466720581055, -0.13620156049728394], [-7.66719913482666, 3.175858736038208, -1.6486923694610596], [5.022003173828125, -4.951065540313721, -1.6797120571136475], [4.285221576690674, -4.214282989501953, -0.41007059812545776], [3.3874499797821045, -3.316511869430542, 0.6318969130516052], [2.363189935684204, -2.2922518253326416, 1.4061485528945923], [1.2518022060394287, -1.180864691734314, 1.8829305171966553], [5.888368129730225, -3.895397901535034, -1.6797120571136475], [5.0220046043396, -3.316511392593384, -0.41007059812545776], [3.966336488723755, -2.6111369132995605, 0.6318969130516052], [2.7619335651397705, -1.8063806295394897, 1.4061485528945923], [1.455080509185791, -0.9331697225570679, 1.8829305171966553], [-1.6190414428710938, 1.0818074941635132, 2.595390796661377], [-3.1758596897125244, 2.1220407485961914, 2.0274128913879395], [-4.610633373260498, 3.080725908279419, 1.1050671339035034], [-5.868222713470459, 3.9210212230682373, -0.13620156049728394], [-6.900296211242676, 4.610631465911865, -1.6486923694610596], [6.532135009765625, -2.6909947395324707, -1.6797120571136475], [5.569483280181885, -2.292250633239746, -0.41007059812545776], [4.396488666534424, -1.8063803911209106, 0.6318969130516052], [3.058227300643921, -1.2520548105239868, 1.4061485528945923], [1.6061298847198486, -0.6505764126777649, 1.8829305171966553], [6.928564548492432, -1.3841410875320435, -1.6797120571136475], [5.906618595123291, -1.1808631420135498, -0.41007059812545776], [4.661374568939209, -0.9331687092781067, 0.6318969130516052], [3.2406840324401855, -0.6505759358406067, 1.4061485528945923], [1.6991455554962158, -0.3439449667930603, 1.8829305171966553], [-1.3768815994262695, 1.3768799304962158, 2.595390796661377], [-2.7008464336395264, 2.7008447647094727, 2.0274128913879395], [-3.921020984649658, 3.921020269393921, 1.1050671339035034], [-4.990512847900391, 4.990513324737549, -0.13620156049728394], [-5.868218898773193, 5.868221282958984, -1.6486923694610596], [7.062422752380371, -0.02505880780518055, -1.6797120571136475], [6.020456314086914, -0.025058681145310402, -0.41007059812545776], [4.750815391540527, -0.025058617815375328, 0.6318969130516052], [3.3022923469543457, -0.025058681145310402, 1.4061485528945923], [1.730553150177002, -0.025058744475245476, 1.8829305171966553], [6.928565502166748, 1.3340235948562622, -1.6797120571136475], [5.906619548797607, 1.1307460069656372, -0.41007059812545776], [4.661374568939209, 0.8830513954162598, 0.6318969130516052], [3.2406842708587646, 0.6004586815834045, 1.4061485528945923], [1.699145793914795, 0.2938274145126343, 1.8829305171966553], [-1.0818090438842773, 1.6190396547317505, 2.595390796661377], [-2.1220414638519287, 3.1758577823638916, 2.0274128913879395], [-3.080726385116577, 4.610631465911865, 1.1050671339035034], [-3.921020746231079, 5.868222713470459, -0.13620156049728394], [-4.610630035400391, 6.90029764175415, -1.6486923694610596], [6.532135963439941, 2.6408774852752686, -1.6797120571136475], [5.569484233856201, 2.242133617401123, -0.41007059812545776], [4.396488666534424, 1.7562637329101562, 0.6318969130516052], [3.058228015899658, 1.2019377946853638, 1.4061485528945923], [1.6061301231384277, 0.6004590392112732, 1.8829305171966553], [5.888369560241699, 3.8452811241149902, -1.6797120571136475], [5.022005081176758, 3.266394853591919, -0.41007059812545776], [3.966336965560913, 2.5610201358795166, 0.6318969130516052], [2.761934280395508, 1.756264328956604, 1.4061485528945923], [1.4550809860229492, 0.883052408695221, 1.8829305171966553], [-0.7451632618904114, 1.798980474472046, 2.595390796661377], [-1.4616879224777222, 3.5288240909576416, 2.0274128913879395], [-2.1220412254333496, 5.12306022644043, 1.1050671339035034], [-2.700845956802368, 6.520419597625732, -0.13620156049728394], [-3.175856351852417, 7.66719913482666, -1.6486923694610596], [5.022005081176758, 4.900949478149414, -1.6797120571136475], [4.28522253036499, 4.164167404174805, -0.41007059812545776], [3.387450695037842, 3.266395330429077, 0.6318969130516052], [2.363190174102783, 2.2421352863311768, 1.4061485528945923], [1.2518031597137451, 1.130747675895691, 1.8829305171966553], [3.966336488723755, 5.767314434051514, -1.6797120571136475], [3.3874504566192627, 4.900949954986572, -0.41007059812545776], [2.6820755004882812, 3.8452818393707275, 0.6318969130516052], [1.877319574356079, 2.6408793926239014, 1.4061485528945923], [1.00410795211792, 1.3340257406234741, 1.8829305171966553], [-0.3798812925815582, 1.9097874164581299, 2.595390796661377], [-0.7451621294021606, 3.7461795806884766, 2.0274128913879395], [-1.0818071365356445, 5.438611030578613, 1.1050671339035034], [-1.3768792152404785, 6.922039985656738, -0.13620156049728394], [-1.6190369129180908, 8.139453887939453, -1.6486923694610596], [2.7619330883026123, 6.411081314086914, -1.6797120571136475], [2.363189697265625, 5.448428630828857, -0.41007059812545776], [1.877319097518921, 4.275433540344238, 0.6318969130516052], [1.3229930400848389, 2.9371728897094727, 1.4061485528945923], [0.7215145826339722, 1.4850751161575317, 1.8829305171966553], [1.4550793170928955, 6.807512283325195, -1.6797120571136475], [1.2518014907836914, 5.78556489944458, -0.41007059812545776], [1.0041069984436035, 4.540319442749023, 0.6318969130516052], [0.721514105796814, 3.1196296215057373, 1.4061485528945923], [0.41488292813301086, 1.5780906677246094, 1.8829305171966553], [0.0959966778755188, 6.941370487213135, -1.6797120571136475], [0.0959966778755188, 5.899401664733887, -0.41007059812545776], [0.0959966778755188, 4.629761219024658, 0.6318969130516052], [0.0959966778755188, 3.1812374591827393, 1.4061485528945923], [0.0959966778755188, 1.609498143196106, 1.8829305171966553]]
shape="rectangle"
def coord(shape):
    if shape=="cone":
        k=volume_cone
    if shape=="rectangle":
        k=Q
    return(k)
    
def mini(liste):
    min_x,min_y,min_z=1000,1000,1000
    for i in liste:
        if i[0]<min_x:
            min_x=i[0]
        if i[1]<min_y:
            min_y=i[1]
        if i[2]<min_z:
            min_z=i[2]
    return([min_x,min_y,min_z])
    
def maxi(liste):
    max_x,max_y,max_z=-1000,-1000,-1000
    for i in liste:
        if i[0]>max_x:
            max_x=i[0]
        if i[1]>max_y:
            max_y=i[1]
        if i[2]>max_z:
            max_z=i[2]
    return([max_x,max_y,max_z])
Xmin,Xmax=mini(coord(shape))[0],maxi(coord(shape))[0]
Ymin,Ymax=mini(coord(shape))[1],maxi(coord(shape))[1]
Zmin,Zmax=mini(coord(shape))[2],maxi(coord(shape))[2]
coord=coord(shape)
k=voroplusplus.compute_voronoi(coord,[[Xmin-1,Xmax+1],[Ymin-1, Ymax+1], [Zmin-1, Zmax+1]],1024.0,radii=[1.3, 1.4])
vertices=[]
faces=[]
#fig = plt.figure()
#ax = plt.axes(projection='3d')
for i in range(len(k)):
    vertices.append(k[i]['vertices'])
for i in range(len(k)):
    faces.append(k[i]['faces'])
for u in range(len(faces)):
    for i in range(len(faces[u])):
        C=faces[u][i]['vertices']
        couples=[]
        for p in range(len(C)-1):
            couples.append([C[p],C[p+1]])
            couples.append([C[-1],C[0]])
            Z=vertices[u]
#        for p in couples:
#           plt.plot([Z[p[0]][0],Z[p[1]][0]],[Z[p[0]][1],Z[p[1]][1]],[Z[p[0]][2],Z[p[1]][2]],'k-',linewidth=0.2)
            




#ax.scatter([1.0, 4.0],[2.0, 5.5],[3.0,6.0],c='b',alpha=0.9)
#for pts in node:
#    ax.scatter(pts[0],pts[1],pts[2],c='b',alpha=0.9)
