// initial
double rootdi = 0.4;
double LAIi = 0.1;

double tilli = 7000.0;
double wrei = 200.0;
double wrti = 4.0;

// functions and parameters for grass
// parameters
double co2a = 360.;  // atmospheric co2 concentration (ppm)
double kdif = 0.60;  //
double LAIcr = 4.;  // critical LAI
double luemax = 3.0;  // light use efficiency (g dm mj-1 par intercepted)
double sla = 0.0025;  // specific leaf area (m2 g-1 dm)
double cLAI = 0.8;  // LAI after cutting (-)
double nitmax = 3.34;  // maximum nitrogen content (%)
double nitr = 3.34;  // actual nitrogen content (%)
double rdrd = 0.01;  // base death rate (fraction)
double tmbas1 = 3.0;  // base temperature perennial ryegrass (degrees celsius)

// parameters for water relations from lingra for thimothee
double drate = 50.0;  // drainage rate (mm day-1)
double irrigf = 0.0;  // irrigation (0 = no irrigation till 1 = full irrigation)
double rootdm = 0.4;  // maximum root depth (m)
double rrdmax = 0.012;  // maximum root growth (m day-1)
double wcad = 0.005;  // air dry water content (fraction)
double wcwp = 0.12;  // wilting point water content (fraction)
double wcfc = 0.29;  // field capacity water content (fraction)
double wci = 0.29;  // initial water content (fraction)
double wcwet = 0.37;  // minimum water content at water logging (fraction)
double wcst = 0.41;  // saturation water content (fraction)

double pi = 3.1415927;
double rad = pi / 180.0;

// initial available water (mm)
double wai = 1000. * rootdi * wci;
// initial leaf weight is initialized as initial
// leaf area divided by initial specific leaf area, kg ha-1
double wlvgi = LAIi / sla;
// remaining leaf weight after cutting is initialized at remaining
// leaf area after cutting divided by initial specific leaf area, kg ha-1
double cwlvg = cLAI / sla;
// maximum site filling new buds (fsmax) decreases due
// to low nitrogen contents, van loo and schapendonk (1992)
// theoretical maximum tillering size = 0.693
double fsmax = nitr / nitmax * 0.693;
