#include <cmath>

double pi = 3.14159265358979;

/*����ģ�ͳ����������ʵ��ֵ��WGS84��*/
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
* DegToRad�����Ƕ�ת��Ϊ���ȡ�
*
*/
inline double DegToRad(double deg){
	return (deg / 180.0 * pi);
}
/*
* UTMCentralMeridian��ȷ������UTM��������������ߡ�
*	Inputs:		zone - ָ��UTM���������ֵ����Χ[1,60]��

*	Returns:	����UTM��������������ߣ���λΪ���ȣ����UTM�������������Χ[1,60]����Ϊ�㡣
				���������ߵķ�Χ���൱��[-177��+177]�Ļ��ȡ�
*/
double UTMCentralMeridian(int zone){
	return DegToRad(-183.0 + (zone * 6.0));
}

/*
* ArcLengthOfMeridian
*	���룺phi�����γ�ȣ���λΪ���ȡ�
*	ȫ��ģ�
*	m_a==����ģ�ͳ��ᡣ
*	sm_b==����ģ�Ͷ��ᡣ
*���أ��㵽�����������룬��λΪ�ס�
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

	������ī����ͶӰ�еľ�γ�ȶ�ת��Ϊx��y���ꡣ
	ע�⣡����ī���У�Transverse Mercator)��UTM��ͬ,֮����Ҫ����������ӡ�

	Inputs:
		phi - ���γ�ȣ���λΪ���ȡ�
	    lambda - ��ľ��ȣ��Ի���Ϊ��λ��
	    lambda0 - Ҫʹ�õ����������ߵľ��ȣ��Ի���Ϊ��λ��
	Outputs:
	    xy - x y-����������x��y�����2Ԫ�����顣
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
	/* ������ķ���ʽ��Ԥ�ȼ���l��ϵ�������������˾Ϳ��Զ�ȡ����ͱ���ı��ʽ
	-- l**1 �� l**2  ���ϵ��1.0 */
	l3coef = 1.0 - t2 + nu2;
	l4coef = 5.0 - t2 + 9 * nu2 + 4.0 * (nu2 * nu2);
	l5coef = 5.0 - 18.0 * t2 + (t2 * t2) + 14.0 * nu2 - 58.0 * t2 * nu2;
	l6coef = 61.0 - 58.0 * t2 + (t2 * t2) + 270.0 * nu2 - 330.0 * t2 * nu2;
	l7coef = 61.0 - 479.0 * t2 + 179.0 * (t2 * t2) - (t2 * t2 * t2);
	l8coef = 1385.0 - 3111.0 * t2 + 543.0 * (t2 * t2) - (t2 * t2 * t2);
	/*���㶫�ࣨx�� */
	xy.x = N * cos(phi) * l + (N / 6.0 * pow(cos(phi), 3.0) * l3coef * pow(l, 3.0))
		+ (N / 120.0 * pow(cos(phi), 5.0) * l5coef * pow(l, 5.0))
		+ (N / 5040.0 * pow(cos(phi), 7.0) * l7coef * pow(l, 7.0));
	/* ���㱱�� (y) */
	xy.y = ArcLengthOfMeridian(phi)
		+ (t / 2.0 * N * pow(cos(phi), 2.0) * pow(l, 2.0))
		+ (t / 24.0 * N * pow(cos(phi), 4.0) * l4coef * pow(l, 4.0))
		+ (t / 720.0 * N * pow(cos(phi), 6.0) * l6coef * pow(l, 6.0))
		+ (t / 40320.0 * N * pow(cos(phi), 8.0) * l8coef * pow(l, 8.0));
}
/*
* LatLonToUTMXY
*	��ͨ�ú���ī����ͶӰ�еľ�γ�ȶ�ת��Ϊx��y���ꡣ
* Inputs:
*   lat - ���γ�ȣ��Ի���Ϊ��λ��
*   lon - ��ľ��ȣ��Ի���Ϊ��λ��
*   zone - ���ڼ���x��yֵ��UTM�����������С��1�����60���򽫸���lon��ֵȷ���ʵ�������
* Outputs:
*   xy - �洢UTM x��yֵ��2Ԫ�����顣
*/

void LatLonToUTMXY(double lat, double lon, int zone, UTMCoor& xy)
{
	MapLatLonToXY(lat, lon, UTMCentralMeridian(zone), xy);

	/*����UTMϵͳ�Ķ���ͱ��ࡣ */
	xy.x = xy.x * UTMScaleFactor + 500000.0;
	xy.y = xy.y * UTMScaleFactor;
	if (xy.y < 0.0)
		xy.y += 10000000.0;
}
/*
	����һ��WGS84��ȥһ��UTM
	��ȷ����Բ����pi = 4 * atan(1.0);
*/
UTMCoor WGS_TO_UTM(WGS84Corr _tempWGS)
{
	pi = 4 * atan(1.0);
	UTMCoor tempUTM;
	// �������
	int nTemp = (int)((_tempWGS.log + 180.0) / 6) + 1;
	LatLonToUTMXY(DegToRad(_tempWGS.lat), DegToRad(_tempWGS.log), nTemp, tempUTM);
	return tempUTM;
}
