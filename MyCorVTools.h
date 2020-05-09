#include <cmath>

double pi = 3.14159265358979;

/*椭球模型常数（这里的实际值是WGS84）*/
double sm_a = 6378137.0;	double sm_b = 6356752.314;
double sm_EccSquared = 6.69437999013e-03;
double UTMScaleFactor = 0.9996;
typedef struct tagUTMCorr{
	double x;
	double y;
}UTMCoor;

typedef struct tagWGS84Corr{
	double lat;
	double log;
}WGS84Corr;
/*
* DegToRad：将角度转换为弧度。
*
*/
inline double DegToRad(double deg){
	return (deg / 180.0 * pi);
}
/*
* UTMCentralMeridian：确定给定UTM区域的中央子午线。
*	Inputs:		zone - 指定UTM区域的整数值，范围[1,60]。

*	Returns:	给定UTM区域的中心子午线，单位为弧度，如果UTM区域参数超出范围[1,60]，则为零。
				中央子午线的范围是相当于[-177，+177]的弧度。
*/
double UTMCentralMeridian(int zone){
	return DegToRad(-183.0 + (zone * 6.0));
}

/*
* ArcLengthOfMeridian
*	输入：phi—点的纬度，单位为弧度。
*	全球的：
*	m_a==椭球模型长轴。
*	sm_b==椭球模型短轴。
*返回：点到赤道的椭球距离，单位为米。
*/

double ArcLengthOfMeridian(double phi)
{
	double alpha, beta, gamma, delta, epsilon, n;
	double result;
	n = (sm_a - sm_b) / (sm_a + sm_b);
	alpha = ((sm_a + sm_b) / 2.0) * (1.0 + (pow(n, 2.0) / 4.0) + (pow(n, 4.0) / 64.0));
	beta = (-3.0 * n / 2.0) + (9.0 * pow(n, 3.0) / 16.0) + (-3.0 * pow(n, 5.0) / 32.0);
	gamma = (15.0 * pow(n, 2.0) / 16.0) + (-15.0 * pow(n, 4.0) / 32.0);
	delta = (-35.0 * pow(n, 3.0) / 48.0) + (105.0 * pow(n, 5.0) / 256.0);
	epsilon = (315.0 * pow(n, 4.0) / 512.0);
	result = alpha * (phi + (beta * sin(2.0 * phi)) + (gamma * sin(4.0 * phi)) 
		+ (delta * sin(6.0 * phi)) + (epsilon * sin(8.0 * phi)));
	return result;
}
/*
* MapLatLonToXY

	将横向墨卡托投影中的经纬度对转换为x和y坐标。
	注意！横轴墨卡托（Transverse Mercator)与UTM不同,之间需要换算比例因子。

	Inputs:
		phi - 点的纬度，单位为弧度。
	    lambda - 点的经度，以弧度为单位。
	    lambda0 - 要使用的中央子午线的经度，以弧度为单位。
	Outputs:
	    xy - x y-包含计算点的x和y坐标的2元素数组。
*/
void MapLatLonToXY(double phi, double lambda, double lambda0, UTMCoor& xy)
{
	double N, nu2, ep2, t, t2, l;
	double l3coef, l4coef, l5coef, l6coef, l7coef, l8coef;
	double tmp;
	ep2 = (pow(sm_a, 2.0) - pow(sm_b, 2.0)) / pow(sm_b, 2.0);
	nu2 = ep2 * pow(cos(phi), 2.0);
	N = pow(sm_a, 2.0) / (sm_b * sqrt(1 + nu2));
	t = tan(phi);
	t2 = t * t;
	tmp = (t2 * t2 * t2) - pow(t, 6.0);
	l = lambda - lambda0;
	/* 在下面的方程式中预先计算l的系数，这样正常人就可以读取东距和北距的表达式
	-- l**1 和 l**2  相关系数1.0 */
	l3coef = 1.0 - t2 + nu2;
	l4coef = 5.0 - t2 + 9 * nu2 + 4.0 * (nu2 * nu2);
	l5coef = 5.0 - 18.0 * t2 + (t2 * t2) + 14.0 * nu2 - 58.0 * t2 * nu2;
	l6coef = 61.0 - 58.0 * t2 + (t2 * t2) + 270.0 * nu2 - 330.0 * t2 * nu2;
	l7coef = 61.0 - 479.0 * t2 + 179.0 * (t2 * t2) - (t2 * t2 * t2);
	l8coef = 1385.0 - 3111.0 * t2 + 543.0 * (t2 * t2) - (t2 * t2 * t2);
	/*计算东距（x） */
	xy.x = N * cos(phi) * l + (N / 6.0 * pow(cos(phi), 3.0) * l3coef * pow(l, 3.0))
		+ (N / 120.0 * pow(cos(phi), 5.0) * l5coef * pow(l, 5.0))
		+ (N / 5040.0 * pow(cos(phi), 7.0) * l7coef * pow(l, 7.0));
	/* 计算北距 (y) */
	xy.y = ArcLengthOfMeridian(phi)
		+ (t / 2.0 * N * pow(cos(phi), 2.0) * pow(l, 2.0))
		+ (t / 24.0 * N * pow(cos(phi), 4.0) * l4coef * pow(l, 4.0))
		+ (t / 720.0 * N * pow(cos(phi), 6.0) * l6coef * pow(l, 6.0))
		+ (t / 40320.0 * N * pow(cos(phi), 8.0) * l8coef * pow(l, 8.0));
}
/*
* LatLonToUTMXY
*	将通用横轴墨卡托投影中的经纬度对转换为x和y坐标。
* Inputs:
*   lat - 点的纬度，以弧度为单位。
*   lon - 点的经度，以弧度为单位。
*   zone - 用于计算x和y值的UTM区域。如果区域小于1或大于60，则将根据lon的值确定适当的区域。
* Outputs:
*   xy - 存储UTM x和y值的2元素数组。
*/

void LatLonToUTMXY(double lat, double lon, int zone, UTMCoor& xy)
{
	MapLatLonToXY(lat, lon, UTMCentralMeridian(zone), xy);

	/*调整UTM系统的东距和北距。 */
	xy.x = xy.x * UTMScaleFactor + 500000.0;
	xy.y = xy.y * UTMScaleFactor;
	if (xy.y < 0.0)
		xy.y += 10000000.0;
}
/*
	进来一个WGS84回去一个UTM
	精确计算圆周率pi = 4 * atan(1.0);
*/
UTMCoor WGS_TO_UTM(WGS84Corr _tempWGS)
{
	pi = 4 * atan(1.0);
	UTMCoor tempUTM;
	// 计算分区
	int nTemp = (int)((_tempWGS.log + 180.0) / 6) + 1;
	LatLonToUTMXY(DegToRad(_tempWGS.lat), DegToRad(_tempWGS.log), nTemp, tempUTM);
	return tempUTM;
}
