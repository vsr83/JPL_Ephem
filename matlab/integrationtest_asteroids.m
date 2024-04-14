% This routine performs numerical integration of the planets in the solar
% system and the Moon 100 years forward and compares the resulting relative
% positions w.r.t. Sun to outputs obtained from the JPL Horizons system

% REFERENCES: 
%  [1] Urban, Seidelmann - Explanatory Supplement to the Astronomical
%  Almanac, 3rd edition, University Science Books, 2013.
%  [2] Horizons System, available at 
%  https://ssd.jpl.nasa.gov/horizons/app.html#/
%  [3] Park, Folkner, Williams, Boggs - The JPL Planetary and Lunar 
%  Ephemerides DE440 and DE441, The Astronomical Journal, vol 161, No 3,
%  2021.

% Initial state of the Moon libration angles at JT_epoch 
% (phi, theta, psi, phi1, theta1, psi1) [1]:
LIB_initial = [
    0.00512995970515812456; 
    0.00004524704499022800;
    0.38239065587686011507; 
    -0.00000223092763198743;
    1.29414222411027863099;
    0.22994485870136698411
];

LIB_initial = [
    0.00512830031411853500;
    0.00004573724185991433;
    0.38239278420173690000;
   -0.00000218986174567295;
    1.29416700274878300000;
    0.22994486018992250000
];

% From DE421 obtained using the python package jplephem.
LIB_initial = [
    5.128132058714363e-03;
    1.165507165777481e-04;
    3.823932005230066e-01;
    1.461912823858170e-05;
    1.294168056057082e+00;
    2.298367282420818e-01
];

% Indices of the bodies.
ind_sun     = 1;
ind_mercury = 2;
ind_venus   = 3;
ind_earth   = 4;
ind_moon    = 5;
ind_mars    = 6;
ind_jupiter = 7;
ind_saturn  = 8;
ind_uranus  = 9;
ind_neptune = 10;
ind_pluto   = 11;
ind_ceres   = 12;
ind_pallas  = 13;
ind_juno    = 14;
ind_vesta   = 15;

% The intial conditions are obtained from JPL Horizons with origin at the
% body center for all bodies except the Sun. This is done since the actual 
% heliocentric initial condition in sources like [1] contain contributions 
% to the Solar System Barycenter from smaller non-planetary objects.

% Initial condition for the position (au).
R_initial = [
    4.501771046200178E-03,  7.667597654825399E-04,  2.662206318863734E-04; ...
    3.572602077668869E-01, -9.154904799747744E-02, -8.598103172768061E-02; ...
    6.082494317526214E-01, -3.491324458403040E-01, -1.955443448720758E-01; ...
    1.160247250568800E-01, -9.265813150908608E-01, -4.017930739420348E-01; ...
    1.152165477160266E-01, -9.285759450811601E-01, -4.028803366239009E-01; ...
   -1.146885853705725E-01, -1.328366526279570E+00, -6.061551990541251E-01;...
   -5.384209277643075E+00, -8.312483870146020E-01, -2.250951187017292E-01; ...
    7.889888161590055E+00,  4.595710989682622E+00,  1.558429793513401E+00; ...
   -1.826990605559244E+01, -1.162723764514111E+00, -2.503714950236171E-01; ...
   -1.605954001605754E+01, -2.394295935145348E+01, -9.400423444456409E+00; ...
   -3.048781548042960E+01, -8.731761132527535E-01,  8.911305385348918E+00; ...
];

% Initial condition for the velocity (au/d).
V_initial = [
   -3.517482063742068E-07,  5.177625639825057E-06,  2.229090837927562E-06; ...
    3.367845709043980E-03,  2.488934292987171E-02,  1.294407129214571E-02; ...
    1.095242018625395E-02,  1.561250662984909E-02,  6.328876605810963E-03; ...
    1.680431652771057E-02,  1.745166161927448E-03,  7.570134278727792E-04; ...
    1.740540134429601E-02,  1.577720694761970E-03,  6.714512881328027E-04; ...
    1.448200480836478E-02,  2.372854532106865E-04, -2.837498250145917E-04; ...
    1.092364404049427E-03, -6.523294106045839E-03, -2.823012134541729E-03; ...
   -3.217204775856843E-03,  4.330632712157481E-03,  1.926417212037553E-03; ...
    2.215425014444429E-04, -3.767652400674729E-03, -1.653244046240149E-03; ...
    2.643121823195712E-03, -1.503490013808834E-03, -6.812710872439278E-04; ...
    3.225591308955338E-04, -3.148753752998151E-03, -1.080178675229524E-03; ...
];

% Astronomical unit (km).
au = 149597870.691;

% Gauss' (gravitational) constant from Table 8.4 in [1].
k = 0.01720209895;

% Earth-Moon mass ratio from Table 8.4 in [1].
em_ratio = 81.30056;

% Standard gravitational parameters from Table 8.4 in [1].
mu = k*k ./ [
    1;
    6023600;
    408523.71;
    328900.5614 * (em_ratio + 1) / em_ratio;
    328900.5614 * (em_ratio + 1);
    3098708;
    1047.3486;
    3497.898;
    22902.98;
    19412.24;
    135200000
]';
%mu(5) =  1.093189443939597e-11
mu(5) = 1.093189450705846e-11;

% From Table 2 in [3].
GM_sun = 132712440041.279419; 
GM_moon= 4902.800118; 
GM_mer = 22031.868551;
GM_ven = 324858.592000; 
GM_ear = 398600.435507;
GM_mar = 42828.375816;
GM_jup = 126712764.100000;
GM_sat = 37940584.841800;
GM_ura = 5794556.400000;
GM_nep = 6836527.100580;
GM_plu = 975.500000;
data_cer        = [62.178, 1.438681557460599E+00, -2.204373520149093E+00, -1.326397869637290E+00,  8.465407071263739E-03,  4.684247756333535E-03,  4.661571593681031E-04];
data_pal        = [13.402, 2.038305484098803E-01, -3.209619169577908E+00,  6.238442521509175E-01,  8.534313988786690E-03, -8.606631066922886E-04, -3.929012607500672E-04];
data_jun        = [1.536, 4.612065292684774E-01, -3.006098944839013E+00, -5.801638036056415E-01,  8.395459067804244E-03,  3.111907833334782E-03,  2.730569008690409E-04];
data_ves        = [17.630, 1.823724941127288E-01,  2.386628127871159E+00,  9.245962963046050E-01, -1.017449623665891E-02,  4.148089092437803E-05,  1.344159277699022E-03];
data_astraea    = [0.159,  2.489298570830910E+00,  1.036393137894717E+00,  2.105623757685725E-01, -5.569106885919010E-03,  7.959737014966118E-03,  3.113961547875325E-03];
data_hebe       = [0.605, 1.339048882519657E+00,  1.442776395814341E+00,  7.927384133525137E-02, -8.775986458159779E-03,  9.426816142936958E-03,  3.535716492239485E-03];
data_iris       = [0.796,  1.892474916134514E+00, -8.484149489370393E-01, -1.571595055065998E-01,  2.786951843966373E-03,  1.131405771472516E-02,  4.975133976616641E-03];
data_flora      = [0.236, -2.119654625312452E+00,  8.084684871277124E-01,  5.333989530870324E-01, -5.818108171130158E-03, -8.811938467706650E-03, -2.835327827441147E-03];
data_metis      = [0.567, -2.424658221530984E+00, -1.253249183418855E-01,  1.859664617477523E-01, -1.166917334659410E-03, -9.845348567586456E-03, -4.559667597322920E-03];
data_hygiea     = [5.364, 2.444258647674053E+00,  2.180590539529486E+00,  1.162854728224503E+00, -5.924502821959189E-03,  5.979689084454833E-03,  2.286440081659974E-03];
data_parthenope = [0.356, -1.231934293166375E+00, -1.941583796431914E+00, -6.486518991650728E-01,  1.014672768156956E-02, -4.310582288320016E-03, -2.342319177719710E-03];
data_egeria     = [0.412, 1.110469088279295E+00, -1.956883565696291E+00, -1.669729882005807E+00,  8.850485918683454E-03,  4.168218701637546E-03,  7.934872925573253E-04];
data_irene      = [0.348, 2.968959645276704E+00,  1.795346737407267E-01, -4.414935256140036E-01,  1.144315321315373E-04,  8.313750673216725E-03,  3.664818354145694E-03];
data_eunomia    = [1.638, -1.438396739682710E+00,  2.001287880286886E+00,  7.672577048909556E-01, -9.735668316806900E-03, -2.981532281300084E-03, -3.694079119333819E-03];
data_psyche     = [2.233, 1.459068533049712E+00, -2.194287373198110E+00, -8.720453549392979E-01,  8.193868988804894E-03,  6.328895956481376E-03,  2.167068152769472E-03];
data_melpomene  = [0.267, -2.742133438940771E+00, -1.276596954534766E-02,  2.398207182934969E-01, -9.963631419254214E-04, -8.953549070722846E-03, -2.237582740412442E-03];
data_fortuna    = [0.463, -2.421847616489512E+00, -1.337430069896114E+00, -5.813969593951555E-01,  4.769320757597536E-03, -7.520732892406689E-03, -2.981884103782758E-03];
data_massalia   = [0.291, -4.469603126635175E-01,  1.855381632768110E+00,  7.780938113421028E-01, -1.251141897166910E-02, -2.534679168776488E-03, -1.141861519707461E-03];
data_lutetia    = [0.139, -4.015026523099796E-01, -2.034736973884179E+00, -8.788573464594445E-01,  1.185954656281440E-02, -1.389992741904971E-04, -7.501838526901043E-04];
data_kalliope   = [0.491, -1.360170232371027E+00, -2.672544296611262E+00, -1.131055313875956E+00,  8.015424090569231E-03, -2.573716161448452E-03, -3.469704531746388E-03];
data_thalia     = [0.129, -1.381182601642768E+00, -2.177230361390210E+00, -8.695469511570262E-01,  7.230264265385984E-03, -5.803213941084448E-03, -4.345035287245586E-03];
data_themis     = [0.403, -1.986326050392570E+00,  1.713081154086989E+00,  7.821777449327153E-01, -7.376544672674655E-03, -7.515391745844364E-03, -3.291592199653334E-03];
data_phocaea    = [0.04, 1.862187673883519E+00, -3.208255256577335E-01,  3.702988878393032E-01,  4.047642354248749E-03,  1.278308592821547E-02,  2.077624549125949E-03];
data_euterpe    = [0.084, -1.260200091759456E+00,  1.418638355224821E+00,  6.489235110334824E-01, -1.081946735425074E-02, -6.731985559609274E-03, -2.575091934900123E-03];
data_bellona    = [0.165, -1.473796820647640E+00,  1.739643324953512E+00,  6.339342179921148E-01, -9.328148846873765E-03, -7.427776303084492E-03, -1.181232128200433E-03];
data_amphitrite = [0.906, -1.469894660268028E-01, -2.362796449275732E+00, -1.340108597116477E+00,  1.008072938451320E-02, -1.922012826080449E-04, -3.957486632603208E-05];
data_urania     = [0.095, -2.606331278449177E+00, -8.031295455197802E-02, -1.194432032888055E-01, -2.716899610319542E-04, -9.151554484547008E-03, -4.227027990320838E-03];
data_euphrosyne = [1.139, -2.344572996538680E+00, -2.202247498588243E+00, -1.499371362567744E+00,  4.929224721470947E-03, -3.894955730352929E-03, -5.777738441181767E-03];
data_daphne     = [0.527, 2.186523624983148E+00,  2.716000195235929E+00,  3.520837863677591E-01, -6.362485913491159E-03,  4.588918690365595E-03,  6.584608776379437E-04];
data_isis       = [0.092, -2.311903896998047E+00,  1.564712483063351E+00,  1.081664679936894E+00, -5.375431543101007E-03, -6.572471474871808E-03, -2.070446989349552E-03];
data_eugenia    = [0.397, -1.286743152948652E+00, -2.064332495462089E+00, -5.815340035179459E-01,  9.624177736943669E-03, -5.447398703614503E-03, -2.359150611235791E-03];
data_nemausa    = [0.144, 2.325725756797201E+00,  9.196308962897707E-01,  1.953459613840103E-01, -4.226339912074973E-03,  9.366246966910217E-03,  2.292475497801344E-03];
data_europa     = [1.354, 1.630502237398202E+00, -2.807392157235787E+00, -1.120111926442581E+00,  7.545905637029689E-03,  4.385623092881690E-03,  6.810559497906901E-04];
data_echo       = [0.021, -2.053757982326704E+00,  5.426195834277842E-01,  1.676053992932709E-01, -5.053759869677038E-03, -1.063488982530823E-02, -3.929518071059859E-03];
data_ausonia    = [0.102, -7.115092271198002E-01, -1.727048748100933E+00, -9.802454407764494E-01,  1.195787236381979E-02, -3.468683811216697E-03, -1.404345165092448E-03];
data_cybele     = [0.694, -2.818819956361019E+00, -1.452782888550498E+00, -4.595465688387552E-01,  5.491205270996683E-03, -7.709769060825189E-03, -2.980911273797440E-03];
data_hesperia   = [0.414, -2.731373448121144E+00, -1.739837131575064E-01, -9.151050728072989E-02, -8.358533004629483E-04, -1.041030935110611E-02, -2.792059302440223E-03];
data_diana      = [0.085, -2.230308906317092E+00, -6.120048256998054E-01, -5.432742940197512E-01,  1.875814559555662E-03, -9.912068098173037E-03, -5.874198924803510E-03];
data_aurora     = [0.414, 1.217207276707023E+00,  2.230212231502064E+00,  1.349204007350788E+00, -9.480914555800900E-03,  3.964169771750117E-03,  2.537799560839387E-03];
data_klotho     = [0.089, -1.909757694987701E+00, -2.678799413270145E+00, -4.520678570587235E-01,  6.291705348125361E-03, -5.039336057956301E-03, -1.538833211441256E-03];
data_ianthe     = [0.055, -2.189380125607045E+00,  2.261161952296765E-01,  1.125269908378213E-01, -2.121674904675346E-03, -9.612548538097386E-03, -7.843087288576463E-03];
data_artemis    = [0.088, -2.150696937406092E+00,  5.418777569041657E-01, -9.913892540157650E-02, -9.340066599988822E-04, -1.186500701241051E-02, -5.102777706698829E-04];
data_ate        = [0.116, -2.065233357643945E+00, -1.395196953626660E+00, -8.511725674554633E-01,  5.618125369204140E-03, -8.131533275794873E-03, -3.600360342257897E-03];
data_hertha     = [0.078, -1.852065892146149E+00, -1.414536880012437E+00, -7.017648232135492E-01,  8.706739887401756E-03, -6.122907894956379E-03, -2.836577493546121E-03];
data_juewa      = [0.188, -2.286169350555527E+00,  1.297951997394274E-01,  1.131118789989825E-01, -1.144931405506583E-03, -1.012917809606640E-02, -6.915212906933010E-03];
data_adeona     = [0.151, 1.101693980964055E+00,  2.046278711333827E+00,  7.335328474510399E-01, -1.047771080417877E-02,  2.702894459461957E-03,  3.883525437203883E-03];
data_lamberta   = [0.105, 1.552360068902740E+00, -1.667018282670181E+00, -1.224506188816949E+00,  1.002692435997029E-02,  4.083988661579011E-03,  1.832028925928610E-03];
data_nausikaa   = [0.107, -2.180257804198296E+00, -1.605282865387759E+00, -1.009730216974740E+00,  6.837111457416405E-03, -5.228655976649039E-03, -2.756824645284942E-03];
data_prokne     = [0.182, 1.495439348517641E+00, -1.277354059314173E+00, -2.998777774810942E-01,  8.680665905284190E-03,  1.043790305312667E-02,  1.341682115827145E-04];
data_kleopatra  = [0.299, -2.623737010950830E+00, -2.141208964541772E+00, -8.484606236444645E-01,  5.369325403114971E-03, -5.867955532138014E-03, -5.850381024164272E-04];
data_athamantis = [0.126, -2.280009064961948E+00, -8.804621513607092E-01, -6.457520952852419E-01,  4.331695284875950E-03, -9.223173324095878E-03, -2.470258241427453E-03];
data_bamberga   = [0.661, 1.398759502667752E+00, -1.287476168320054E+00, -6.690981722049137E-01,  7.164361562381234E-03,  9.219960292753789E-03,  6.857862930338518E-03];
data_devosa     = [0.033, 2.057439975917922E+00, -1.322163020004281E+00, -7.809757236714602E-01,  5.188028454852560E-03,  7.591643630714724E-03,  4.670779725911250E-03];
data_desiderata = [0.114, -1.436277775979161E+00,  2.243809450657797E+00,  2.055782163205734E+00, -6.473878619730446E-03, -4.308989152437247E-03, -1.132191205623079E-03];
data_eleonora   = [0.327, -4.005999335179206E-01, -2.845495704101578E+00, -3.632228628478518E-01,  9.411660735247279E-03, -2.169502767815410E-03, -2.270976648186976E-03];
data_palma      = [0.355, -2.501243697853252E+00, -1.473357710199123E+00, -2.186350250404082E+00,  4.414112690330698E-03, -5.694497887780149E-03, -4.127569372791904E-03];
data_thia       = [0.092, -1.680815360167212E+00, -7.209597941389104E-01, -6.362643130177040E-01,  6.053480195657973E-03, -1.200452216953581E-02, -3.193084806523710E-03];
data_aspasia    = [0.216, 2.627942559142085E+00,  3.285668199803020E-02,  4.993790865859914E-01, -1.211476343770462E-04,  9.800913393411879E-03,  3.214451230789752E-03];
data_aurelia    = [0.102, -8.322059812298748E-01,  2.925043492697467E+00,  1.070612510793479E+00, -7.845381414317704E-03, -2.434066420406861E-03, -1.376246240291387E-03];
data_patientia  = [0.61, 1.259011084632854E+00,  2.454040194218057E+00,  6.904822353910874E-01, -9.063538223399895E-03,  3.448055766345343E-03,  4.183421795838043E-03];
data_kreusa     = [0.164, -2.067224908799407E+00, -1.883391717418183E+00, -3.931589343775042E-01,  5.745444200230644E-03, -7.721755199184345E-03, -4.781672319841375E-03];
data_davida     = [1.638, -2.160191690204074E+00,  1.486363084137172E+00,  1.096959372997934E+00, -7.390695339994901E-03, -7.779023127591244E-03, -4.855410222679363E-04];
data_herculina  = [0.886, -2.931622799810343E-01, -2.481690401594876E+00, -7.259337890912891E-01,  1.015343129054321E-02, -2.073913322188581E-03, -3.634396471848115E-03];
data_peraga     = [0.044, 1.567515380989014E+00, -1.756839186460132E+00, -7.301456430553257E-01,  7.168202539106875E-03,  7.126475338220530E-03,  3.646324489422235E-03];
data_zelinda    = [0.09, 2.454317195892013E+00, -1.156415988994922E+00,  3.118162946236253E-01,  1.774158702407912E-03,  7.979676512769986E-03,  4.592146350030268E-03];
data_interamnia = [2.464, 2.462834949672756E+00, -1.150472583413863E-01,  7.842197898051008E-01, -1.309572315126532E-03,  1.034898731363664E-02,  4.812523167671943E-03];
data_winchester = [0.196, -8.349019593402328E-01,  2.059697397518895E+00,  6.132931334064005E-01, -1.213473829480265E-02, -2.004138913456460E-03,  2.575813568314367E-03];

asteroid_data = [
    data_cer;
    data_pal;
    data_jun;
    data_ves;
    data_astraea;
    data_hebe;
    data_iris;
    data_flora;
    data_metis;
    data_hygiea;
    data_parthenope;
    data_egeria;
    data_irene;
    data_eunomia;
    data_psyche;
    data_melpomene;
    data_fortuna;
    data_massalia;
    data_lutetia;
    data_kalliope;
    data_thalia;
    data_themis;
    data_phocaea;
    data_euterpe;
    data_bellona;
    data_amphitrite;
    data_urania;
    data_euphrosyne;
    data_daphne;
    data_isis;
    data_eugenia;
    data_nemausa;
    data_europa;
    data_echo;
    data_ausonia;
    data_cybele;
    data_hesperia;
    data_diana;
    data_aurora;
    data_klotho;
    data_ianthe;
    data_artemis;
    data_ate;
    data_hertha;
    data_juewa;
    data_adeona;
    data_lamberta;
    data_nausikaa;
    data_prokne;
    data_kleopatra;
    data_athamantis;
    data_bamberga;
    data_devosa;
    data_desiderata;
    data_eleonora;
    data_palma;
    data_thia;
    data_aspasia;
    data_aurelia;
    data_patientia;
    data_kreusa;
    data_davida;
    data_herculina;
    data_peraga;
    data_zelinda;
    data_interamnia;
    data_winchester
];

mu(2) = k * k * GM_mer / GM_sun;
mu(3) = k * k * GM_ven / GM_sun;
mu(4) = k * k * GM_ear / GM_sun;
mu(5) = k * k * GM_moon / GM_sun;
mu(6) = k * k * GM_mar / GM_sun;
mu(7) = k * k * GM_jup / GM_sun;
mu(8) = k * k * GM_sat / GM_sun;
mu(9) = k * k * GM_ura / GM_sun;
mu(10)= k * k * GM_nep / GM_sun;
mu(11)= k * k * GM_plu / GM_sun;
%mu(12)= k * k * GM_cer / GM_sun;
%mu(13)= k * k * GM_pal / GM_sun;
%mu(14)= k * k * GM_jun / GM_sun;
%mu(15)= k * k * GM_ves / GM_sun;

num_asteroids = length(asteroid_data);
num_targets = 11 + num_asteroids;

for ind_asteroid = 1:num_asteroids
    mu(11 + ind_asteroid) = k * k * asteroid_data(ind_asteroid, 1) / GM_sun;
    R_initial(11 + ind_asteroid, :) = asteroid_data(ind_asteroid, 2:4);
    V_initial(11 + ind_asteroid, :) = asteroid_data(ind_asteroid, 5:7);
end

% Load data for checking the accuracy of the results.
jpl_reference;

% Make coordinates consistent with the Sun expressed with respect to SSB by 
% adding it to all OSVs.
for ind_target = 2:num_targets    
    R_initial(ind_target, :) = R_initial(ind_target, :) + R_initial(1, :);
    V_initial(ind_target, :) = V_initial(ind_target, :) + V_initial(1, :);
end

% OSVs related to the DoFs during the integration.
R = R_initial;
V = V_initial;

% Compute relativistic barycenter.
[r_bary, v_bary] = barycenter(R, V, mu, true)

% Offset bodies w.r.t. so that the relativistic barycenter is approximately
% at the origin.
for ind_target = 1:num_targets
    R(ind_target, :) = R(ind_target, :) - r_bary;
    V(ind_target, :) = V(ind_target, :) - v_bary;
end

reference = {
    [0, 0, 0, 0, 0];
    mercury_exp;
    venus_exp;
    earth_exp;
    moon_exp;
    mars_exp;
    jupiter_exp;
    saturn_exp;
    uranus_exp;
    neptune_exp;
    pluto_exp
};
names = ["Sun", "Mercury", "Venus", "Earth", "Moon", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune", "Pluto"];
names_short = ["S", "M", "V", "E", "Mo", "Ma", "J", "S", "U", "N", "P"];


% Matrix with position errors in (km) for all targets for each timestep.
ERR = [];
% Time step size in days.
JT_timestep = 0.01;
% Start date 1969-Jun-28 00:00:00.0000
JT_epoch = 2440400.50;
% Compute number of timesteps.
num_timesteps = floor(365.25 * 100 / JT_timestep + 1);
% Current Julian time.
JT = JT_epoch;
% Time after epoch.
t = 0;


[R, V] = moon_frombary(R, V, mu, 4, 5);

% We first have 6 DoFs for the libration angles and their time derivatives
% and then additional 6 DoFs (x, y, z, x1, y1, z1) for each body.
y = [LIB_initial; osv_to_dof(R, V)];
[R, V] = moon_tobary(R, V, mu, 4, 5);

y_start2 = y;
R_start2 = R;
V_start2 = V;


ERR = [];
LIB_out = [];
JT = t + JT_epoch;

for timestep = 1:num_timesteps
    % Prepare output print to display during the integration.
    err_step = zeros(1, num_targets);
    err_str = "";
    filled = false;

    for ind_target = 1:num_targets
        if ind_target > ind_pluto
            continue
        end

        ref_target = reference{ind_target};
        results = find(abs(ref_target(:, 1) - JT) < 1/86400);
        if length(results) > 0
            year         = ref_target(results, 2);
            if ~filled
                err_str = sprintf("%d - %.10f ", year, JT);
            end

            r_target     = R(ind_target, 1:3) - R(1, 1:3);
            r_target_exp = ref_target(results, 3:5);
            err_target = norm(r_target - r_target_exp) * au;

            err_step(1, ind_target) = err_target;

            err_str_target = sprintf("%s %.2f ", names_short(ind_target), err_target);
            err_str = strcat(err_str, err_str_target);
            filled = true;
        end        
    end
    if filled
        disp(err_str);
        ERR = [ERR; year, err_step];
        LIB_out = [LIB_out; y(1:6)'];
    end

    if timestep < 5
        %disp(y(1:6)');
    end

    % Perform the actual integration.
    [y, t] = runge4(@func_lib, t, y, JT_timestep, mu);

    % Convert DoFs back to OSVs.
    [R, V] = dof_to_osv(y(7:end));
    [R, V] = moon_tobary(R, V, mu, 4, 5);

    JT = t + JT_epoch;

    if filled
        % Visualize the error.
        figure(1);
        clf
        
        year        = ERR(:, 1);
        err_mercury = ERR(:, 3);
        err_venus   = ERR(:, 4);
        err_earth   = ERR(:, 5);
        err_moon    = ERR(:, 6);
        err_mars    = ERR(:, 7);
        err_jupiter = ERR(:, 8);
        err_saturn  = ERR(:, 9);
        err_uranus  = ERR(:, 10);
        err_neptune = ERR(:, 11);
        err_pluto   = ERR(:, 12);
        
        endyear = 1969 + length(ERR) - 1;
        semilogy(year, err_mercury, 'r',   'LineWidth', 2)
        hold on;
        semilogy(year, err_venus,   'g',   'LineWidth', 2)
        semilogy(year, err_earth,   'b',   'LineWidth', 2)
        semilogy(year, err_moon,    'bo',  'LineWidth', 2)
        semilogy(year, err_mars,    'm',   'LineWidth', 2)
        semilogy(year, err_jupiter, 'c',   'LineWidth', 2)
        semilogy(year, err_saturn,  'y',   'LineWidth', 2)
        semilogy(year, err_uranus,  'k',   'LineWidth', 2)
        semilogy(year, err_neptune, 'k--', 'LineWidth', 2)
        semilogy(year, err_pluto,   'k:',  'LineWidth', 2)
        
        legend('Mercury', 'Venus', 'Earth', 'Moon', 'Mars', 'Jupiter', 'Saturn', 'Uranus', ...
            'Neptune', 'Pluto', 'FontSize', 18);
        xlabel('Year', 'FontSize', 18);
        ylabel('Error in Position (km)', 'FontSize', 18);
        xlim([1969, endyear])
        ylim([0.01 1000])
        grid on
        title(sprintf('Numerical Integration Error %d-%d', min(year), max(year)), 'FontSize', 18);
    end
end