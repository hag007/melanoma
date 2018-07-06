import scipy.special
import matplotlib.pyplot as plt
from matplotlib import style
style.use("ggplot")
import scipy
from scipy.stats import hypergeom
from statsmodels.sandbox.stats.multicomp import fdrcorrection0
import logging
sh = logging.StreamHandler()
logger = logging.getLogger("log")
logger.addHandler(sh)
from constants import *
from infra import *

# pvals = [0.9813, 0.4420, 0.2890, 0.9791, 0.9991, 0.7856, 0.0314, 0.6535, 0.1499, 0.2626, 0.0437]
#counts
# pvals = [2.72E-01, 9.73E-01, 3.27E-06, 3.79E-02, 2.18E-01, 9.09E-01]
# fpkm
pvals = [0.0004120246766123503, 0.0005946290792310789, 0.0008977000862616192, 0.0010693695725064786, 0.001500821476386621, 0.0015692037169598905, 0.001698230092281935, 0.0030522409215740995, 0.003272082329470733, 0.0035574156822292904, 0.00469868420974227, 0.0049018846399946575, 0.004963479809638217, 0.005815441581120635, 0.00610234663778063, 0.0062289753088665386, 0.006719985553708915, 0.0075594305280683866, 0.007635807965803652, 0.007865363258428338, 0.008232027317389168, 0.008304045640007427, 0.008464993845824567, 0.00890851706394925, 0.009265297883270813, 0.009370874788704017, 0.009881943485523275, 0.009952890493257446, 0.01047122247816675, 0.011170745399583562, 0.012351051173722023, 0.012440555199137913, 0.012525115632071007, 0.012675263129740249, 0.013348962368628023, 0.013708564744434389, 0.014256736634076443, 0.014656779262308565, 0.01583224320730394, 0.016095415607341533, 0.01649298649188556, 0.016587253500627563, 0.01672164012154446, 0.01748595237098747, 0.017739325394859343, 0.019084469332838708, 0.019807840063882507, 0.020355787093707837, 0.020425565881554003, 0.020997818378226403, 0.02102521692381843, 0.02112850655047881, 0.02121638790764658, 0.021720571691580374, 0.021898213711872977, 0.02219763240221762, 0.02263905663291071, 0.022719259596244, 0.02320616646513626, 0.023606941391503367, 0.023882697473384493, 0.024028682232829238, 0.024617536151323532, 0.024893432790539884, 0.02539301162309144, 0.025511970733047257, 0.02566695505801727, 0.025800020934005515, 0.026236130247321302, 0.026384578558772555, 0.026843626954602934, 0.027059741778198802, 0.02747259900698002, 0.027679326251889107, 0.027997049199902396, 0.028149799038818618, 0.029057402105270427, 0.029377076839829784, 0.029739777593583018, 0.029780245315818916, 0.03001249521136688, 0.030090389893966243, 0.03182456630959991, 0.031850978059240205, 0.03202050897584454, 0.033874687273682036, 0.03399774985465395, 0.0340866287482396, 0.03512878902897565, 0.03559941563446986, 0.03594744831842575, 0.03640234291873921, 0.03770396715737052, 0.03898498072902997, 0.04028439173606591, 0.04112779375610101, 0.04168810764608577, 0.04169182506114214, 0.0423428161952366, 0.04260297370696926, 0.04272943044261874, 0.043460567910168266, 0.043592304875088105, 0.04367498784800605, 0.0449858406348382, 0.045200538678651994, 0.045401234126930146, 0.04635331060544971, 0.04660914116982807, 0.04750847888256012, 0.048052012271036545, 0.048111575183611935, 0.04822631599480463, 0.0484439405899416, 0.04869737446750875, 0.048778472605343265, 0.050223939263078776, 0.05022676492577818, 0.050569693793956406, 0.050676799511957804, 0.05068577422148586, 0.05074074113093025, 0.050848487227668944, 0.05185785871641503, 0.0532134005993931, 0.053340860024512535, 0.05349063627775593, 0.054006519140617, 0.05417211929029225, 0.05610786938531135, 0.05663723151551348, 0.05714208267277795, 0.05718476016264617, 0.0585352212812221, 0.058730897090082515, 0.059706116132087615, 0.059955639412658676, 0.05997643072385187, 0.060332488731715345, 0.060484911540903075, 0.062486338601391356, 0.06344046074653262, 0.06425092760830166, 0.06539349009839052, 0.06540420268966907, 0.0673297891976032, 0.06733274476119096, 0.06779976219511176, 0.07066143340153161, 0.07087343590263684, 0.07178962658309042, 0.07222833361171213, 0.07243365517458243, 0.07299275607082364, 0.0734291633354059, 0.07392208668251363, 0.07423028917802227, 0.07489467570597148, 0.07526430822432771, 0.07536897251707161, 0.07775125748808107, 0.07824431929018755, 0.07999404747761862, 0.08016307061665325, 0.08083457928877406, 0.08106265410638867, 0.0811259033378135, 0.08158433185262766, 0.08162380766749064, 0.08223081131133993, 0.08277624735578062, 0.08284902227473812, 0.08310479541762222, 0.08326260796371006, 0.0832761840138041, 0.0841596327907256, 0.08424383102178948, 0.08442083651972586, 0.0850041785566745, 0.08502550094508304, 0.08645925373339138, 0.0869174113140034, 0.08929198190589224, 0.0899808363013486, 0.09046992761275402, 0.09053030383573715, 0.09066421339939774, 0.09147971549582147, 0.0916188964259523, 0.09302250801771093, 0.09373790263578889, 0.09437409574413354, 0.09538231974487732, 0.09577145893045534, 0.09628836434672648, 0.09634919847302874, 0.09640000533415657, 0.09667577614389726, 0.09695300235200942, 0.09729952498790996, 0.10110314257001986, 0.1014285078744363, 0.10169900548892634, 0.10346794897257505, 0.10386351495146122, 0.10412771543819195, 0.10414610519815883, 0.10484338665839842, 0.10491631583432094, 0.10718000695955789, 0.10720194199441058, 0.10737329825294133, 0.10739067216064829, 0.10740927622433445, 0.10757702389031018, 0.10766012543045285, 0.10912409611579874, 0.10937434930427986, 0.10991984361221721, 0.11247895193983808, 0.11304025909710783, 0.11383471675184349, 0.11429695787273873, 0.11477712577985921, 0.11480239597999436, 0.11520059628765862, 0.11697108340903087, 0.11783681885402178, 0.11846025667958916, 0.11847379529163939, 0.11931151322745186, 0.11984122827820497, 0.12011658074595997, 0.12049117168674553, 0.12191567655279892, 0.12200221317198823, 0.1221127738320307, 0.1227855544909234, 0.12294339861159337, 0.12311213474678263, 0.12392304600978643, 0.12468511758346539, 0.12469602355445275, 0.12516721593678912, 0.12523387464905722, 0.12570917851175487, 0.12572817679716128, 0.12602265482380262, 0.12671420778788253, 0.12749847355576427, 0.1281690769915898, 0.1287510531196883, 0.12876561648166743, 0.12983851148468512, 0.13036057897663644, 0.1310004932539362, 0.13399902296354882, 0.1345462389402872, 0.1347107421567832, 0.13478879685196726, 0.135159445530543, 0.13529163644599987, 0.13540992814339506, 0.13558035372401384, 0.13573225420684054, 0.1358490578719882, 0.13653378688554127, 0.1366088992072568, 0.1369052539824776, 0.1370696695509149, 0.13741196141128506, 0.13827757002336194, 0.13875768504222522, 0.13922856622829483, 0.1403188651540967, 0.14201835076006528, 0.14225889441797496, 0.1426054331182033, 0.1431821004646717, 0.14435200475270438, 0.1443629280410505, 0.14554291203507924, 0.14581978808401322, 0.14618267180868608, 0.1470260915844841, 0.1471201540689715, 0.14761121674860406, 0.14782316694852193, 0.14831103717657185, 0.14959626482446078, 0.15007704318559822, 0.15074308266103448, 0.15216264598637663, 0.15216942385580548, 0.15243057125998416, 0.15277219323442545, 0.15372369144169837, 0.1540021088835455, 0.1544672106714145, 0.15452722820075332, 0.15482221874758495, 0.15554960976765675, 0.15579941870526845, 0.15635927092578727, 0.15708844661829566, 0.15728889750687028, 0.1576963051972444, 0.1578435285297433, 0.1578462700485138, 0.1616736777699469, 0.1617381313826789, 0.1618525419587485, 0.1618682271932198, 0.16213620546196247, 0.1621538365425209, 0.16285457386634147, 0.1641692948291692, 0.16432880535673308, 0.1655738140538054, 0.16568092681590696, 0.1662617100847197, 0.1675113651743019, 0.16892629155688876, 0.1691602162584146, 0.16924729243291545, 0.1705470685988156, 0.1705828843232663, 0.1708210136293091, 0.17203843253331275, 0.17225934231752646, 0.1735485825696077, 0.1742443269351293, 0.1742918768591701, 0.17434046976788956, 0.17612053639041264, 0.1766594802968452, 0.17677462419852052, 0.17771707928654112, 0.17823476747421915, 0.1795854111638036, 0.17972885016188292, 0.17973697850203507, 0.17974488503304625, 0.17977467699473668, 0.1799163937134942, 0.1799451402617436, 0.1799841248523388, 0.18006657663835365, 0.18008403382913707, 0.18070787083348336, 0.18103269729913032, 0.18139880664350425, 0.1818672634206461, 0.18215672580835415, 0.18271437943955787, 0.18295878764871829, 0.18318975299812923, 0.18393143587204372, 0.1842964829058447, 0.18449648887783848, 0.18462053160457692, 0.18472813702758548, 0.1850278651891392, 0.18504383448711773, 0.18522914360682782, 0.1859369847074451, 0.1873285621792068, 0.1876848196259015, 0.1886857027702155, 0.18877761387571979, 0.18883321338479245, 0.18909599563311388, 0.18944724727847026, 0.19148599817796247, 0.19177678380099467, 0.19197339434567404, 0.19404902071311234, 0.19524100171239595, 0.19540334495939465, 0.1955645271213662, 0.19655400830689848, 0.19694708726747995, 0.19722701976978477, 0.19770497918047208, 0.19777578967660478, 0.19830577726894869, 0.19896565735733301, 0.19921834027605648, 0.2003675875872667, 0.20058731613104247, 0.200599329797355, 0.20224047402484288, 0.20292437857041476, 0.20336706697979565, 0.2053424969985521, 0.2060836413937335, 0.2073124140477499, 0.20757076442828423, 0.20933192179367532, 0.20957787792293217, 0.20971065470645311, 0.2113968433408401, 0.21201104648559935, 0.21233962920401744, 0.21246901757183145, 0.21355475440093707, 0.2155794414672801, 0.21672556560518783, 0.2171138972262269, 0.2172734342091678, 0.21923939551808638, 0.21926269539997495, 0.21938985389953894, 0.2208557093349283, 0.22203100685417834, 0.22226558497993693, 0.2223625260618001, 0.22250183639906598, 0.22334198908004996, 0.2233594277362333, 0.22382824629448247, 0.22396172102633113, 0.22462228025633083, 0.225164528168952, 0.22596312518487655, 0.22638466189655176, 0.2269726581751689, 0.22726491610187544, 0.2286860215416203, 0.22935922840491332, 0.23072892122282546, 0.23279994225103726, 0.2330347716783619, 0.23567637376447784, 0.23616578362567134, 0.2373448514290285, 0.2373842531202513, 0.2385528031780473, 0.2390168967914771, 0.23992547317434532, 0.24236432863618745, 0.24262242370700382, 0.24274616072807337, 0.24605387462480857, 0.24715292431762448, 0.24771779930233678, 0.24890935795690705, 0.2492879100243945, 0.2503699925511872, 0.2505145775842837, 0.2508054427641709, 0.25135131974854863, 0.25172688839064067, 0.2518273051248936, 0.2523592420964288, 0.25390630605558356, 0.25540643299299665, 0.25551542830236196, 0.2557741163834009, 0.256223634476177, 0.2572926387174618, 0.257651357043461, 0.25923804727973554, 0.2596846767347784, 0.259792663049787, 0.26111001615158613, 0.2634135201071853, 0.2634661484078157, 0.2635116464564697, 0.26375751907964695, 0.26575494828491913, 0.2679899071255801, 0.26878499204308504, 0.26919790530766635, 0.2694083026674365, 0.2701894390643099, 0.2710137272648973, 0.2767824164499864, 0.27758977362173093, 0.2801665280532071, 0.28027422016817644, 0.280734701651292, 0.2832580925567021, 0.2851186003890687, 0.2879393011992594, 0.2884705341864864, 0.2902461115487734, 0.2925031694707754, 0.2930787663259609, 0.29365766297810414, 0.29365766297810414, 0.29365766297810414, 0.29365766297810414, 0.29365766297810414, 0.29365766297810414, 0.29365766297810414, 0.29365766297810414, 0.2936576629781042, 0.2936576629781042, 0.2936576629781042, 0.2936576629781042, 0.2936576629781042, 0.2936576629781042, 0.2936576629781042, 0.2936576629781042, 0.2936576629781042, 0.2936576629781042, 0.2936576629781042, 0.2936576629781042, 0.2936576629781042, 0.2936576629781042, 0.2936576629781042, 0.2936576629781042, 0.2936576629781042, 0.2936576629781042, 0.2936576629781042, 0.2936576629781042, 0.2936576629781042, 0.2936576629781042, 0.2936576629781042, 0.2936576629781042, 0.2936576629781042, 0.2936576629781042, 0.2936576629781042, 0.2936576629781042, 0.2936576629781042, 0.2936576629781042, 0.29365766297810436, 0.29365766297810436, 0.29365766297810436, 0.29365766297810436, 0.29365766297810436, 0.29365766297810436, 0.29365766297810436, 0.29365766297810436, 0.29365766297810436, 0.29365766297810436, 0.29365766297810436, 0.29365766297810436, 0.2959007381699327, 0.2963402092188547, 0.2969025534357258, 0.29866587295647473, 0.30063938892745395, 0.30499021307358753, 0.3058333085676798, 0.306597330024437, 0.30666791277928596, 0.307063162364396, 0.30767542134232473, 0.30825901476492096, 0.3090057763703929, 0.31061571175369657, 0.31093396207811835, 0.3117600163219818, 0.31274529776447674, 0.31295802839808484, 0.31296054793004485, 0.314089836991997, 0.3147797535370396, 0.3164742703888551, 0.3176551948543012, 0.31841847523927275, 0.3188678646250245, 0.32153334584548227, 0.32154463814337, 0.3216212774456549, 0.3220538643922396, 0.32315502180931643, 0.32321856440999686, 0.32431651270739603, 0.3244919694381896, 0.325253239515116, 0.32573420180586654, 0.3268257056656384, 0.32684036857551935, 0.32711386152028854, 0.32862798371085145, 0.32881250336941925, 0.32930329815176995, 0.32930329815176995, 0.3305721200096422, 0.33178626185556204, 0.33233478096154423, 0.33291967718031035, 0.33330566789918525, 0.3349994485253942, 0.3351377277514437, 0.33570932245242435, 0.3368517284000694, 0.3373014989495636, 0.33776961654655346, 0.3385725319315579, 0.3392120086382069, 0.3399727065842558, 0.34059502269982667, 0.34224441330615163, 0.34292902389588487, 0.3443659329337042, 0.344685396901266, 0.34474607806599855, 0.34541256270737497, 0.34772643381755264, 0.34780500086052624, 0.34781062403242646, 0.34781062403242646, 0.34781062403242646, 0.34781062403242646, 0.34781062403242646, 0.34781062403242646, 0.34781062403242646, 0.34781062403242646, 0.34781062403242646, 0.34781062403242646, 0.34781062403242646, 0.34781062403242646, 0.34781062403242646, 0.34781062403242646, 0.34781062403242646, 0.34781062403242646, 0.34781062403242646, 0.34781062403242646, 0.34781062403242646, 0.34781062403242646, 0.34781062403242646, 0.34781062403242646, 0.34781062403242646, 0.34781062403242646, 0.34781062403242646, 0.34781062403242646, 0.34781062403242646, 0.34781062403242646, 0.34781062403242646, 0.34781062403242646, 0.34781062403242646, 0.34781062403242646, 0.34781062403242646, 0.3478106240324267, 0.3478106240324267, 0.3478106240324267, 0.3478106240324267, 0.3478106240324267, 0.3478106240324267, 0.3478106240324267, 0.3478106240324267, 0.3478106240324267, 0.3478106240324267, 0.3478106240324267, 0.3478106240324267, 0.34850915787036263, 0.3488738901704008, 0.3513675762289846, 0.3523171705776421, 0.3552304551847002, 0.35556808897611947, 0.3570141003620052, 0.35707886809357814, 0.35736750594621813, 0.357816708123763, 0.3603383029641366, 0.3620191737847842, 0.36378990492143803, 0.364320668282753, 0.3650533691253619, 0.36589991868394034, 0.36615818979817616, 0.3663021750284442, 0.36716291978099547, 0.36814376077818645, 0.3687526615404281, 0.36888366079254475, 0.3699661691498314, 0.3700067165538484, 0.3702809439065581, 0.372135935006918, 0.3729567540582095, 0.3729951132874756, 0.3735009532361384, 0.3735518779117689, 0.373663865525554, 0.37375272563740647, 0.37445381050643267, 0.37506489049963776, 0.37578100254015145, 0.37678090628134475, 0.38097566386537485, 0.38187082025744823, 0.38240248477349537, 0.38338951145002653, 0.38344471875694885, 0.3837106071577576, 0.3839818949602273, 0.38466860905402656, 0.38473472780357587, 0.3852445889427675, 0.38596220067087206, 0.3863644566373946, 0.38800091633353107, 0.3889743331519595, 0.3893007968699357, 0.3904327066322796, 0.39129787376919944, 0.39174205356678915, 0.39222958049945744, 0.3930831788124417, 0.3939942449035789, 0.39500979446491824, 0.3966821501057064, 0.39702117952059535, 0.39714871621870773, 0.3975171896646382, 0.39768959185666486, 0.39809832717023985, 0.3988101554507252, 0.39893457585694325, 0.39974198461812094, 0.39978571844452493, 0.40124868303529615, 0.40138706663327817, 0.4017083228790542, 0.40182642265202917, 0.4056158697273554, 0.4056291026816058, 0.4068000053570281, 0.4070965744582017, 0.4084692519684294, 0.40898073190676354, 0.4094365347273824, 0.4096331516718418, 0.41205950822439263, 0.4125550267141914, 0.4127854245872473, 0.412914150350695, 0.4132518729314717, 0.4146839294268291, 0.41695257458544155, 0.41787328529711354, 0.4188293666953681, 0.41926598593915054, 0.41933823940407733, 0.41979048217868664, 0.4198014535821437, 0.4228472455699207, 0.4245379924869881, 0.4247694024532259, 0.4249392106311205, 0.42514795217965506, 0.42595378983169585, 0.42654574197447814, 0.4268259333841605, 0.42727020283225325, 0.4278104464558654, 0.42791735028978517, 0.4308104226342202, 0.4309796388718313, 0.4327907325106679, 0.4335447608834643, 0.43381165376557385, 0.4342290839972325, 0.43472326470701106, 0.43474561095963493, 0.4373518224768097, 0.43833061767408, 0.43857066850648874, 0.4394926242075211, 0.43987396542474944, 0.44042476256313934, 0.4407256777297931, 0.4422254182926998, 0.44271293817006896, 0.44290935009515653, 0.4434190027267465, 0.44383021985887683, 0.44547653868460335, 0.4455706432417703, 0.4459915315035977, 0.4471410513001457, 0.44809018182712035, 0.45037223424914563, 0.45141613444461803, 0.4520366844376206, 0.4528985081499023, 0.45293717872419936, 0.45449152127450565, 0.4547205015051603, 0.4549083285130441, 0.45506630694851324, 0.4560408441477385, 0.45689789716404305, 0.4572226589606385, 0.4573315264676575, 0.45802438207697704, 0.4585050959365289, 0.4594962252185746, 0.4604043044164815, 0.4608112060607241, 0.4623013218585236, 0.46356108137278385, 0.4635662886395977, 0.46374656479856646, 0.46534311300706177, 0.46573688471004404, 0.4660582084501438, 0.46679486772965884, 0.4685365292003435, 0.4698877932950767, 0.46994512687383305, 0.470622379161859, 0.4714915167984135, 0.47168978667507055, 0.47348613326846356, 0.4739966519523817, 0.47470630820142645, 0.4761354638532118, 0.4771337573586333, 0.4776155206841234, 0.4789028642873969, 0.4789430993266174, 0.4795872163251108, 0.47959716952913456, 0.4796802200301906, 0.48074583406980953, 0.48087989442082346, 0.4809930105400002, 0.4810482792609919, 0.48183434675273307, 0.4825500593030593, 0.4838933909383295, 0.48442348940459634, 0.4857176665234685, 0.486299936323522, 0.4867577287259882, 0.4869153924865479, 0.4874396380585998, 0.4875180098120244, 0.4889420398215364, 0.4912349661887131, 0.4913837813948818, 0.494094181824034, 0.4958336066180522, 0.4963569765429996, 0.4967944873952048, 0.49709501846016146, 0.49885332993870923, 0.49931056793478246, 0.49979584240240893, 0.500650470156857, 0.5011957517644097, 0.5025212612476986, 0.5049208829450933, 0.5050821134027177, 0.5059698701988302, 0.507032096854009, 0.5070420990088498, 0.5072376016104425, 0.5083245314524699, 0.5092702666881244, 0.5095507785594324, 0.5097778134554856, 0.5097850184600246, 0.5100194977477839, 0.5100978458736286, 0.5102411764625219, 0.5108630427561991, 0.5110652811496907, 0.5116608998679842, 0.5123577275680957, 0.5145776674649609, 0.5155249801352486, 0.5162866205369422, 0.5168943712984541, 0.5169033424394446, 0.516904233804787, 0.5178124341642789, 0.5181320763975115, 0.518357795073741, 0.5190383813982491, 0.5197119014482148, 0.520393122353813, 0.5229882931032509, 0.5232317805873578, 0.5234028257301514, 0.524162400069838, 0.524344546388023, 0.5246038737925323, 0.5254480767626512, 0.5254907866886096, 0.5264038160491422, 0.5264546694129448, 0.526777150168088, 0.5301244878478111, 0.5302569798174395, 0.5323359809432039, 0.5331867457035818, 0.5339606792197875, 0.5340335528940165, 0.535117410981214, 0.5365679727137522, 0.5366972009242087, 0.5377829176073694, 0.5389802685197036, 0.5401448830856507, 0.5405564662871556, 0.5414568138087481, 0.5419789849766032, 0.5422656162799611, 0.5426333931505649, 0.5449919372246478, 0.5455587636276467, 0.5461600219596119, 0.5466906268905343, 0.5468198974475795, 0.5469452974350977, 0.5472164768736292, 0.5484282136890036, 0.5487996648406241, 0.5514272212932101, 0.5518772307849708, 0.5525221945778142, 0.5528273621733758, 0.5533131571810004, 0.5584207611088711, 0.5603193325803433, 0.5611202949783144, 0.5612288629572582, 0.5614200972706249, 0.5623599522545489, 0.5633684416086719, 0.5650894497510778, 0.5651041808614805, 0.5669277561138582, 0.5672783608092842, 0.5673689004192689, 0.5681963355525549, 0.5685826302941398, 0.5694783659659036, 0.5758006702810325, 0.5779909532539067, 0.5780880623579894, 0.57847867273832, 0.5790245056921206, 0.5797437421346304, 0.5801133154137165, 0.5814063336681816, 0.5815517601516527, 0.5823580951095186, 0.5823957400796491, 0.5875024198679494, 0.5884782062259796, 0.5889364166541435, 0.5893901597564222, 0.5894382200999639, 0.5896625392761528, 0.5899372111922013, 0.5902541148559586, 0.590461091328643, 0.5907215921939251, 0.5914832667919229, 0.5918084578021866, 0.5938219077715496, 0.5947828682312861, 0.5949099869780388, 0.596908308273684, 0.5989916571106747, 0.6006754873874391, 0.6012423371410249, 0.6023416914922981, 0.6031623898340244, 0.6040772932610323, 0.6055401233067628, 0.6061528639021696, 0.6067917486049326, 0.6082433036382661, 0.6087643856003639, 0.6093629356831909, 0.6099159586173689, 0.6109568768532094, 0.6110694796056051, 0.6127438792972625, 0.6131312124205911, 0.6141369171356572, 0.6154757098348209, 0.616303743990954, 0.6165481989369392, 0.6168421786312513, 0.6179947759057098, 0.6186225274686443, 0.6194845934858162, 0.6211228879216963, 0.6242814681955288, 0.625115625573389, 0.6297421191513882, 0.6314650659636429, 0.6320594301867055, 0.6333979542219144, 0.6334415407960141, 0.6348850893941347, 0.6359109015985374, 0.6367151180551509, 0.6391503879481485, 0.6401140735308937, 0.6402485652181373, 0.6425107534899858, 0.6432724107926668, 0.6437302072172495, 0.6443878940665024, 0.6447244341021607, 0.6447904236193513, 0.6451982474380052, 0.6452351722947702, 0.6456757702867304, 0.6461180356982024, 0.647652458866923, 0.6479785762670729, 0.6490014618341478, 0.6504722576488757, 0.6507975527667696, 0.6520339668023108, 0.6532235290345637, 0.6543959967479803, 0.6555467009798268, 0.6582107119134366, 0.6590855939565756, 0.6592600595286808, 0.6646568468862829, 0.6661556494044922, 0.66623594243701, 0.6666062185939012, 0.668613388742896, 0.6687136024379199, 0.6694433124106762, 0.6705425922147399, 0.6710298765605159, 0.6718145827880103, 0.674566666618268, 0.6745872328364746, 0.6746228643360282, 0.6753420204786225, 0.6769247232469823, 0.6780148320093788, 0.6811588422637551, 0.6819057283462613, 0.682311587679002, 0.6828378579093939, 0.6834012598416462, 0.6844373765995941, 0.6847321925961067, 0.6849662691681055, 0.6851203521994298, 0.6860608202124286, 0.687655435868405, 0.6887174703340306, 0.6897615392033672, 0.6915264120844733, 0.6919038172609058, 0.6920776433430949, 0.6923753003614468, 0.6947401217488363, 0.6967116037840745, 0.6974483105524268, 0.698156697404401, 0.6988610198935461, 0.6997113782554064, 0.702341955980188, 0.7036133801344796, 0.7038623921203765, 0.7041369966483724, 0.7047953308649562, 0.7054742470961513, 0.7065380121361708, 0.707043872808089, 0.7084405168508545, 0.7104931500257262, 0.7112350940390693, 0.7116299038191309, 0.7148992718429411, 0.7151247867089006, 0.7154398089697782, 0.7179428393101648, 0.7199741891277293, 0.7200517267510178, 0.7201590508987361, 0.7204130700693305, 0.7212741778614942, 0.7223697229432755, 0.7237879999344463, 0.724175756734754, 0.7246631002837352, 0.7257365958529772, 0.7258169410627755, 0.7263743216994583, 0.7267951298555622, 0.7289849001983948, 0.7290428007525624, 0.7300374591203769, 0.731494719428218, 0.7327536180751537, 0.7333860236093925, 0.7348066612221358, 0.7349737007182828, 0.737474854917687, 0.7376663473429573, 0.7388315219709447, 0.73994805042806, 0.7410230681190182, 0.7411993875104217, 0.7422428496941421, 0.7427971421475998, 0.7440706029547426, 0.7445504946741961, 0.7454181572456362, 0.746201607068063, 0.7478622439522847, 0.7496525584681665, 0.7502718807848624, 0.7523621003478465, 0.7536267476695335, 0.7540140503258063, 0.7560525710442073, 0.7568297842664274, 0.7571607477658229, 0.7575162805464868, 0.7578210910715544, 0.7579090617667885, 0.7581214841900472, 0.758911755381498, 0.7599239617010722, 0.7599332403403224, 0.7609473217222269, 0.7611682070035649, 0.761927208375497, 0.7621862256146379, 0.7636869399115654, 0.7638104132078887, 0.7644993794000984, 0.7645271133956781, 0.7647829995590414, 0.7652157361158574, 0.7661922098703708, 0.7662745077610628, 0.7667698345165437, 0.7682101002177933, 0.7702064573599681, 0.7703393516199394, 0.7707603987081655, 0.7708573126926201, 0.7728439606455793, 0.7751427856218347, 0.7754378770818784, 0.7767937194333351, 0.7773572781281074, 0.7779003704769704, 0.778806288936212, 0.7792769413105434, 0.7801087356728343, 0.7804064080157301, 0.7804566767961179, 0.7826804133490395, 0.7842491352120876, 0.7885872629270112, 0.7887707679882946, 0.7894417480409858, 0.7906153844847934, 0.7917164444529562, 0.7917840225844414, 0.7925252207650783, 0.7925707953096544, 0.7926344833004833, 0.7936386325829821, 0.7941797922541243, 0.7967201619252967, 0.7967533304275879, 0.79853425539737, 0.8004999854062722, 0.8013593936759914, 0.8055250055044596, 0.8079743097842134, 0.8083268120902247, 0.809500630285767, 0.8109805802090564, 0.8133403918639472, 0.8155112206574419, 0.8155450486771627, 0.8157927356695749, 0.8166650787865495, 0.8188995960360487, 0.8210835132683162, 0.8221818455804845, 0.8231650062769558, 0.8231846797254406, 0.8232359521789725, 0.823692206400326, 0.8254487563178703, 0.8261588848688871, 0.8282862922073003, 0.8284094591953024, 0.8309085997815913, 0.8318291128706118, 0.8320865663736697, 0.8340959962118838, 0.8350664929843217, 0.835613362617539, 0.8358250104905112, 0.8358812921792089, 0.8393389238277063, 0.8407044043873906, 0.8409745577841996, 0.8410558455152668, 0.8435619802513705, 0.8444478079058241, 0.8452013833240781, 0.8468626391418483, 0.8471265994500775, 0.8479541461583101, 0.848043451617829, 0.8481249029750316, 0.8483106593501979, 0.8484523026887596, 0.8494861449782806, 0.8499307693022982, 0.8501551450447586, 0.8508134912365084, 0.8521294198309617, 0.8561025694772089, 0.8572163036429156, 0.8595519847342372, 0.8603390360760605, 0.8615262715239931, 0.8618413051433058, 0.8618488509213175, 0.8628547021860182, 0.8632030328418581, 0.8644922427779822, 0.8646299493755685, 0.8657903066915622, 0.8662987544219924, 0.8684887036378881, 0.8702389088037009, 0.8710088534539211, 0.8710566868711825, 0.8729789183343845, 0.873275920957463, 0.8749194404312617, 0.8756165906585259, 0.8760639416674426, 0.8767488414917, 0.8768911274037664, 0.8784703811269431, 0.8790259727809125, 0.8795797319617285, 0.8796782578961193, 0.8816573503087277, 0.8817586521909266, 0.8829422358161433, 0.8851697693555786, 0.8852284318864434, 0.8866371926898183, 0.8874368067907165, 0.8899741583009976, 0.8901019042604498, 0.8901677679131978, 0.8908453941718344, 0.8909178965888769, 0.8910080190404774, 0.8918695497206789, 0.891977523495041, 0.892528520397648, 0.8935364891273168, 0.8948186599980259, 0.8954934627831201, 0.8963148510569918, 0.897704324509669, 0.8987182164361667, 0.9000890523590792, 0.9008525043222804, 0.9015823537002223, 0.9017212984225612, 0.9019620635591601, 0.9021453937950308, 0.904899495973198, 0.9053164607011975, 0.9053298081213349, 0.9055405588831269, 0.9058646873219856, 0.9060629402355971, 0.9108689621574044, 0.9144094556471899, 0.9191166830179998, 0.9196598126169626, 0.9201920130510287, 0.9206199536973456, 0.9211225706316687, 0.9211244766108373, 0.9214256247046901, 0.9215221650855789, 0.9216489642423836, 0.9216630439931817, 0.9224642796347815, 0.9226057291405447, 0.9226279384427399, 0.9234364632647812, 0.9243721116070343, 0.9250736860635741, 0.9253832593487041, 0.9259784492432147, 0.9273450706214443, 0.9285994637329298, 0.9286979282332688, 0.9296605949900099, 0.9299070050058599, 0.9308147570918113, 0.9309879872952942, 0.9310871635001353, 0.9322800885430198, 0.9325861911730063, 0.9332916274357264, 0.934304649747317, 0.9343286473229888, 0.9348251625489691, 0.9351573578467296, 0.9362946913768636, 0.9372539035906137, 0.940548379551214, 0.9413999911561569, 0.942172117059388, 0.9437912158959758, 0.9444769991219414, 0.9453259190811577, 0.9453743421620585, 0.945849065838821, 0.9459387666853101, 0.9490835017923489, 0.9496458162558145, 0.9499250024305024, 0.9501686564171272, 0.9507759826705571, 0.9508077722684964, 0.9511466878786468, 0.9513635905843546, 0.9520457886259241, 0.9522871732357141, 0.9537623648868385, 0.9543828925625899, 0.9550120112016702, 0.956439948262535, 0.9572619042042421, 0.9574498938714305, 0.9583388317598576, 0.9605773999231089, 0.9607178066120434, 0.9614158168586834, 0.9615389770928655, 0.9627035191648666, 0.9629072197827988, 0.9631747294774198, 0.9642953885920126, 0.9659514281830439, 0.9674936415384291, 0.9676923047571331, 0.9684096919642654, 0.9695207030750981, 0.9698080782135451, 0.9736847712324288, 0.9750301927427982, 0.9803391029847579, 0.9804311520752395, 0.9827642940645644, 0.9832199114808234, 0.9851046153202172, 0.9851326634170271, 0.9864613108890081, 0.9867606753255681, 0.9867983853846606, 0.9877214582237583, 0.9884814030398696, 0.9893361128870052, 0.9919789299685842, 0.9934678575841119, 0.9955652258984142, 0.9955933191746463, 0.9956643864189734, 0.9967313476204092, 0.9977487197567388, 0.9989128515772294]
# fpkm-uq
# pvals = [0.04364191582855, 0.0746803276988296, 3.23745624934477E-09, 0.0337253006483806, 0.00202830018481174, 0.0278338429119658]
pvals.sort()
fdr_results = fdrcorrection0(pvals, alpha=0.05, method='indep', is_sorted=True)
true_counter = len([cur for cur in fdr_results[0] if cur == True])
print fdr_results[1]
print "true hypothesis: {}".format(true_counter)
print "total hypothesis: {}".format(np.size(fdr_results[0]))