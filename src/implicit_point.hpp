/****************************************************************************
* Indirect predicates for geometric constructions					        *
*                                                                           *
* Consiglio Nazionale delle Ricerche                                        *
* Istituto di Matematica Applicata e Tecnologie Informatiche                *
* Sezione di Genova                                                         *
* IMATI-GE / CNR                                                            *
*                                                                           *
* Authors: Marco Attene                                                     *
* Copyright(C) 2019: IMATI-GE / CNR                                         *
* All rights reserved.                                                      *
*                                                                           *
* This program is free software; you can redistribute it and/or modify      *
* it under the terms of the GNU Lesser General Public License as published  *
* by the Free Software Foundation; either version 3 of the License, or (at  *
* your option) any later version.                                           *
*                                                                           *
* This program is distributed in the hope that it will be useful, but       *
* WITHOUT ANY WARRANTY; without even the implied warranty of                *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser  *
* General Public License for more details.                                  *
*                                                                           *
* You should have received a copy of the GNU Lesser General Public License  *
* along with this program.  If not, see http://www.gnu.org/licenses/.       *
*                                                                           *
****************************************************************************/

#include "implicit_point.h"
#include "indirect_predicates.h"

int orient2d(double p1x, double p1y, double p2x, double p2y, double p3x, double p3y);
int orient3d(double px, double py, double pz, double qx, double qy, double qz, double rx, double ry, double rz, double sx, double sy, double sz);

template<class PT, class T>
inline bool lambda2d_SSI(
	const PT* l1, const PT* l2, const PT* m1, const PT* m2,
	T& lambda_x, T& lambda_y, T& lambda_det)
{
	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundUP();

	const T t1 = l1[0] * l2[1] - l2[0] * l1[1];
	const T t3 = m1[0] * m2[1] - m2[0] * m1[1];
	const T tx2 = m1[0] * m2[2] - m2[0] * m1[2];
	const T tx4 = l1[0] * l2[2] - l2[0] * l1[2];
	const T ty2 = m1[1] * m2[2] - m2[1] * m1[2];
	const T ty4 = l1[1] * l2[2] - l2[1] * l1[2];
	lambda_x = t1 * tx2 - t3 * tx4;
	lambda_y = t1 * ty2 - t3 * ty4;
	lambda_det = tx4 * ty2 - tx2 * ty4;

	if constexpr (std::is_same<interval_number, T>::value) {
		setFPUModeToRoundNEAR();
		return lambda_det.signIsReliable();
	}
	else return true;
}

template<class PT, class T>
inline bool lambda3d_LPI(const PT* p, const PT* q, const PT* r, const PT* s, const PT* t, T& lambda_x, T& lambda_y, T& lambda_z, T& lambda_d) {
	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundUP();

	const T trz = (t[2] * r[3] - r[2] * t[3]);
	const T sry = (s[1] * r[3] - r[1] * s[3]);
	const T trY = (t[1] * r[3] - r[1] * t[3]);
	const T srx = (s[0] * r[3] - r[0] * s[3]);
	const T trx = (t[0] * r[3] - r[0] * t[3]);
	const T srz = (s[2] * r[3] - r[2] * s[3]);

	const T a2233 = trz * sry - trY * srz;
	const T a2133 = trz * srx - trx * srz;
	const T a2132 = trY * srx - trx * sry;

	const T pxqd = p[0] * q[3];
	const T qxpd = q[0] * p[3];
	const T pyqd = p[1] * q[3];
	const T qypd = q[1] * p[3];
	const T pzqd = p[2] * q[3];
	const T qzpd = q[2] * p[3];
	const T pdqd = p[3] * q[3];

	const T comm = (pxqd - qxpd) * a2233 - (pyqd - qypd) * a2133 + (pzqd - qzpd) * a2132;
	const T conn = (pxqd * r[3] - r[0] * pdqd) * a2233 + (pzqd * r[3] - r[2] * pdqd) * a2132 - (pyqd * r[3] - r[1] * pdqd) * a2133;
	const T corr = comm * r[3];

	lambda_x = corr * pxqd - (pxqd - qxpd) * conn;
	lambda_y = corr * pyqd - (pyqd - qypd) * conn;
	lambda_z = corr * pzqd - (pzqd - qzpd) * conn;
	lambda_d = corr * pdqd;

	if constexpr (std::is_same<interval_number, T>::value) {
		setFPUModeToRoundNEAR();
		return lambda_d.signIsReliable();
	}
	else return true;
}

template<class PT, class T>
inline bool lambda3d_TPI(const PT* v1, const PT* v2, const PT* v3, const PT* w1, const PT* w2, const PT* w3, const PT* u1, const PT* u2, const PT* u3, T& lambda_x, T& lambda_y, T& lambda_z, T& lambda_d) {
	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundUP();

	const T v12y = v1[3] * v2[1] - v1[1] * v2[3];
	const T v23z = v2[3] * v3[2] - v2[2] * v3[3];
	const T v12z = v1[3] * v2[2] - v1[2] * v2[3];
	const T v23y = v2[3] * v3[1] - v2[1] * v3[3];
	const T v23x = v2[3] * v3[0] - v2[0] * v3[3];
	const T v12x = v1[3] * v2[0] - v1[0] * v2[3];
	const T nvx0 = v12y * v23z - v12z * v23y;
	const T nvy0 = v12z * v23x - v12x * v23z;
	const T nvz0 = v12x * v23y - v12y * v23x;
	const T w12y = w1[3] * w2[1] - w1[1] * w2[3];
	const T w23z = w2[3] * w3[2] - w2[2] * w3[3];
	const T w12z = w1[3] * w2[2] - w1[2] * w2[3];
	const T w23y = w2[3] * w3[1] - w2[1] * w3[3];
	const T w23x = w2[3] * w3[0] - w2[0] * w3[3];
	const T w12x = w1[3] * w2[0] - w1[0] * w2[3];
	const T nwx0 = w12y * w23z - w12z * w23y;
	const T nwy0 = w12z * w23x - w12x * w23z;
	const T nwz0 = w12x * w23y - w12y * w23x;
	const T u12y = u1[3] * u2[1] - u1[1] * u2[3];
	const T u23z = u2[3] * u3[2] - u2[2] * u3[3];
	const T u12z = u1[3] * u2[2] - u1[2] * u2[3];
	const T u23y = u2[3] * u3[1] - u2[1] * u3[3];
	const T u23x = u2[3] * u3[0] - u2[0] * u3[3];
	const T u12x = u1[3] * u2[0] - u1[0] * u2[3];
	const T nux0 = u12y * u23z - u12z * u23y;
	const T nuy0 = u12z * u23x - u12x * u23z;
	const T nuz0 = u12x * u23y - u12y * u23x;
	const T nv01 = nvx0 * v1[0] + nvy0 * v1[1] + nvz0 * v1[2];
	const T nw01 = nwx0 * w1[0] + nwy0 * w1[1] + nwz0 * w1[2];
	const T nu01 = nux0 * u1[0] + nuy0 * u1[1] + nuz0 * u1[2];
	const T nvwu = nv01 * w1[3] * u1[3];
	const T nuvw = nu01 * v1[3] * w1[3];
	const T nwvu = nw01 * v1[3] * u1[3];
	lambda_x = nvwu * (nwy0 * nuz0 - nwz0 * nuy0) + nuvw * (nvy0 * nwz0 - nvz0 * nwy0) - nwvu * (nvy0 * nuz0 - nvz0 * nuy0);
	lambda_y = nwvu * (nvx0 * nuz0 - nvz0 * nux0) - nvwu * (nwx0 * nuz0 - nwz0 * nux0) - nuvw * (nvx0 * nwz0 - nvz0 * nwx0);
	lambda_z = nuvw * (nvx0 * nwy0 - nvy0 * nwx0) + nvwu * (nwx0 * nuy0 - nwy0 * nux0) - nwvu * (nvx0 * nuy0 - nvy0 * nux0);
	lambda_d = (nvx0 * (nwy0 * nuz0 - nwz0 * nuy0) + nvz0 * (nwx0 * nuy0 - nwy0 * nux0) - nvy0 * (nwx0 * nuz0 - nwz0 * nux0)) * v1[3] * w1[3] * u1[3];

	if constexpr (std::is_same<interval_number, T>::value) {
		setFPUModeToRoundNEAR();
		return lambda_d.signIsReliable();
	}
	else return true;
}

template<class PT, class UT, class T>
inline bool lambda3d_LNC(const PT* p, const PT* q, const UT t, T& lambda_x, T& lambda_y, T& lambda_z, T& lambda_d) {
	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundUP();

	const T a = (1.0 - t) * q[3];
	const T b = p[3] * t;
	lambda_x = p[0] * a + q[0] * b;
	lambda_y = p[1] * a + q[1] * b;
	lambda_z = p[2] * a + q[2] * b;
	lambda_d = p[3] * q[3];

	if constexpr (std::is_same<interval_number, T>::value) {
		setFPUModeToRoundNEAR();
		return lambda_d.signIsReliable();
	}
	else return true;
}

template<class PT, class UT, class T>
inline bool lambda3d_BPT(const PT* p, const PT* q, const PT* r, const UT u, const UT v, T& lambda_x, T& lambda_y, T& lambda_z, T& lambda_d) {
	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundUP();

	const T rdqd = r[3] * q[3];
	const T pdrdu = p[3] * r[3] * u;
	const T ccc = p[3] * q[3] * (v + u - 1.0);
	const T rdvqd = rdqd * v;
	lambda_x = q[0] * pdrdu - r[0] * ccc + p[0] * rdvqd;
	lambda_y = q[1] * pdrdu - r[1] * ccc + p[1] * rdvqd;
	lambda_z = q[2] * pdrdu - r[2] * ccc + p[2] * rdvqd;
	lambda_d = p[3] * rdqd;

	if constexpr (std::is_same<interval_number, T>::value) {
		setFPUModeToRoundNEAR();
		return lambda_d.signIsReliable();
	}
	else return true;
}

template<class PT, class T>
inline bool lambda3d_TBC(const PT* p, const PT* q, const PT* r, const PT* s, T& lambda_x, T& lambda_y, T& lambda_z, T& lambda_d) {
	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundUP();

	lambda_d = p[3] * q[3] * r[3] * s[3];
	const T dpqr = p[3] * q[3] * r[3];
	const T dpqs = p[3] * q[3] * s[3];
	const T dprs = p[3] * r[3] * s[3];
	const T dqrs = q[3] * r[3] * s[3];
	lambda_x = (dqrs * p[0] + dprs * q[0] + dpqs * r[0] + dpqr * s[0]) * 0.25;
	lambda_y = (dqrs * p[1] + dprs * q[1] + dpqs * r[1] + dpqr * s[1]) * 0.25;
	lambda_z = (dqrs * p[2] + dprs * q[2] + dpqs * r[2] + dpqr * s[2]) * 0.25;

	if constexpr (std::is_same<interval_number, T>::value) {
		setFPUModeToRoundNEAR();
		return lambda_d.signIsReliable();
	}
	else return true;
}

inline int lessThan_EE(double x1, double y1, double z1, double x2, double y2, double z2)
{
	int ret;
	if ((ret = ((x1 > x2) - (x1 < x2)))) return ret;
	if ((ret = ((y1 > y2) - (y1 < y2)))) return ret;
	return ((z1 > z2) - (z1 < z2));
}

inline int lessThan_IE(const genericPoint& p1, double x2, double y2, double z2)
{
	int ret;
	if ((ret = lessThanOnX_IE(p1, x2))) return ret;
	if ((ret = lessThanOnY_IE(p1, y2))) return ret;
	return lessThanOnZ_IE(p1, z2);
}

inline int lessThan_II(const genericPoint& p1, const genericPoint& p2)
{
	int ret;
	if ((ret = lessThanOnX_II(p1, p2))) return ret;
	if ((ret = lessThanOnY_II(p1, p2))) return ret;
	return lessThanOnZ_II(p1, p2);
}

inline int lessThan_EE(const genericPoint& a, const genericPoint& b) { return lessThan_EE(a.toExplicit3D().X(), a.toExplicit3D().Y(), a.toExplicit3D().Z(), b.toExplicit3D().X(), b.toExplicit3D().Y(), b.toExplicit3D().Z()); }
inline int lessThan_IE(const genericPoint& a, const genericPoint& b) { return lessThan_IE(a, b.toExplicit3D().X(), b.toExplicit3D().Y(), b.toExplicit3D().Z()); }

inline int genericPoint::lessThan(const genericPoint& a, const genericPoint& b)
{
	if (a.isExplicit3D() && b.isExplicit3D()) return lessThan_EE(a, b);
	if (!a.isExplicit3D() && b.isExplicit3D()) return lessThan_IE(a, b);
	if (a.isExplicit3D() && !b.isExplicit3D()) return -lessThan_IE(b, a);
	return lessThan_II(a, b);
}

inline int lessThanOnX_IE(const genericPoint& a, const genericPoint& b) { return lessThanOnX_IE(a, b.toExplicit3D().X()); }

inline int genericPoint::lessThanOnX(const genericPoint& a, const genericPoint& b)
{
	if (a.isExplicit3D() && b.isExplicit3D())
	{
		double av = a.toExplicit3D().X(), bv = b.toExplicit3D().X();
		return ((av > bv) - (av < bv));
	}

	if (!a.isExplicit3D() && b.isExplicit3D()) return lessThanOnX_IE(a, b);
	if (a.isExplicit3D() && !b.isExplicit3D()) return -lessThanOnX_IE(b, a);
	return lessThanOnX_II(a, b);
}

inline int lessThanOnY_IE(const genericPoint& a, const genericPoint& b) { return lessThanOnY_IE(a, b.toExplicit3D().Y()); }

inline int genericPoint::lessThanOnY(const genericPoint& a, const genericPoint& b)
{
	if (a.isExplicit3D() && b.isExplicit3D())
	{
		double av = a.toExplicit3D().Y(), bv = b.toExplicit3D().Y();
		return ((av > bv) - (av < bv));
	}
	if (!a.isExplicit3D() && b.isExplicit3D()) return lessThanOnY_IE(a, b);
	if (a.isExplicit3D() && !b.isExplicit3D()) return -lessThanOnY_IE(b, a);
	return lessThanOnY_II(a, b);
}

inline int lessThanOnZ_IE(const genericPoint& a, const genericPoint& b) { return lessThanOnZ_IE(a, b.toExplicit3D().Z()); }

inline int genericPoint::lessThanOnZ(const genericPoint& a, const genericPoint& b)
{
	if (a.isExplicit3D() && b.isExplicit3D())
	{
		double av = a.toExplicit3D().Z(), bv = b.toExplicit3D().Z();
		return ((av > bv) - (av < bv));
	}
	if (!a.isExplicit3D() && b.isExplicit3D()) return lessThanOnZ_IE(a, b);
	if (a.isExplicit3D() && !b.isExplicit3D()) return -lessThanOnZ_IE(b, a);
	return lessThanOnZ_II(a, b);
}

inline int dotproductSign2D_EEE(const genericPoint& a, const genericPoint& b, const genericPoint& c) { return dotProductSign2D(a.toExplicit2D().X(), a.toExplicit2D().Y(), b.toExplicit2D().X(), b.toExplicit2D().Y(), c.toExplicit2D().X(), c.toExplicit2D().Y()); }
inline int dotproductSign2D_IEE(const genericPoint& a, const genericPoint& b, const genericPoint& c) { return dotProductSign2D_IEE(a, b.toExplicit2D().X(), b.toExplicit2D().Y(), c.toExplicit2D().X(), c.toExplicit2D().Y()); }
inline int dotproductSign2D_IIE(const genericPoint& a, const genericPoint& b, const genericPoint& c) { return dotProductSign2D_IIE(a, b, c.toExplicit2D().X(), c.toExplicit2D().Y()); }
inline int dotproductSign2D_III(const genericPoint& a, const genericPoint& b, const genericPoint& c) { return dotProductSign2D_III(a, b, c); }
inline int dotproductSign2D_EEI(const genericPoint& a, const genericPoint& b, const genericPoint& c) { return dotProductSign2D_EEI(c, a.toExplicit2D().X(), a.toExplicit2D().Y(), b.toExplicit2D().X(), b.toExplicit2D().Y()); }
inline int dotproductSign2D_IEI(const genericPoint& a, const genericPoint& b, const genericPoint& c) { return dotProductSign2D_IEI(a, c, b.toExplicit2D().X(), b.toExplicit2D().Y()); }

inline int genericPoint::dotProductSign2D(const genericPoint& a, const genericPoint& b, const genericPoint& c)
{
	if (a.isExplicit2D() && b.isExplicit2D() && c.isExplicit2D()) return dotproductSign2D_EEE(a, b, c);
	if (a.isSSI() && b.isExplicit2D() && c.isExplicit2D()) return dotproductSign2D_IEE(a, b, c);
	if (a.isExplicit2D() && b.isSSI() && c.isExplicit2D()) return dotproductSign2D_IEE(b, a, c);
	if (a.isExplicit2D() && b.isExplicit2D() && c.isSSI()) return dotproductSign2D_EEI(a, b, c);
	if (a.isSSI() && b.isSSI() && c.isExplicit2D()) return dotproductSign2D_IIE(a, b, c);
	if (a.isExplicit2D() && b.isSSI() && c.isSSI())	return dotproductSign2D_IEI(b, a, c);
	if (a.isSSI() && b.isExplicit2D() && c.isSSI())	return dotproductSign2D_IEI(a, b, c);
	return dotproductSign2D_III(a, b, c);
}

inline int dotproductSign3D_EEE(const genericPoint& a, const genericPoint& b, const genericPoint& c) { return dotProductSign3D(a.toExplicit3D().X(), a.toExplicit3D().Y(), a.toExplicit3D().Z(), b.toExplicit3D().X(), b.toExplicit3D().Y(), b.toExplicit3D().Z(), c.toExplicit3D().X(), c.toExplicit3D().Y(), c.toExplicit3D().Z()); }
inline int dotproductSign3D_IEE(const genericPoint& a, const genericPoint& b, const genericPoint& c) { return dotProductSign3D_IEE(a, b.toExplicit3D().X(), b.toExplicit3D().Y(), b.toExplicit3D().Z(), c.toExplicit3D().X(), c.toExplicit3D().Y(), c.toExplicit3D().Z()); }
inline int dotproductSign3D_IIE(const genericPoint& a, const genericPoint& b, const genericPoint& c) { return dotProductSign3D_IIE(a, b, c.toExplicit3D().X(), c.toExplicit3D().Y(), c.toExplicit3D().Z()); }
inline int dotproductSign3D_III(const genericPoint& a, const genericPoint& b, const genericPoint& c) { return dotProductSign3D_III(a, b, c); }
inline int dotproductSign3D_EEI(const genericPoint& a, const genericPoint& b, const genericPoint& c) { return dotProductSign3D_EEI(c, a.toExplicit3D().X(), a.toExplicit3D().Y(), a.toExplicit3D().Z(), b.toExplicit3D().X(), b.toExplicit3D().Y(), b.toExplicit3D().Z()); }
inline int dotproductSign3D_IEI(const genericPoint& a, const genericPoint& b, const genericPoint& c) { return dotProductSign3D_IEI(a, c, b.toExplicit3D().X(), b.toExplicit3D().Y(), b.toExplicit3D().Z()); }

inline int genericPoint::dotProductSign3D(const genericPoint& a, const genericPoint& b, const genericPoint& c)
{
	if (a.isExplicit3D() && b.isExplicit3D() && c.isExplicit3D()) return dotproductSign3D_EEE(a, b, c);
	if (!a.isExplicit3D() && b.isExplicit3D() && c.isExplicit3D()) return dotproductSign3D_IEE(a, b, c);
	if (a.isExplicit3D() && !b.isExplicit3D() && c.isExplicit3D()) return dotproductSign3D_IEE(b, a, c);
	if (a.isExplicit3D() && b.isExplicit3D() && !c.isExplicit3D()) return dotproductSign3D_EEI(a, b, c);
	if (!a.isExplicit3D() && !b.isExplicit3D() && c.isExplicit3D()) return dotproductSign3D_IIE(a, b, c);
	if (a.isExplicit3D() && !b.isExplicit3D() && !c.isExplicit3D())	return dotproductSign3D_IEI(b, a, c);
	if (!a.isExplicit3D() && b.isExplicit3D() && !c.isExplicit3D())	return dotproductSign3D_IEI(a, b, c);
	return dotproductSign3D_III(a, b, c);
}

inline int orient2d_EEE(const genericPoint& a, const genericPoint& b, const genericPoint& c) { return orient2d(a.toExplicit2D().X(), a.toExplicit2D().Y(), b.toExplicit2D().X(), b.toExplicit2D().Y(), c.toExplicit2D().X(), c.toExplicit2D().Y()); }
inline int orient2d_IEE(const genericPoint& a, const genericPoint& b, const genericPoint& c) { return orient2d_indirect_IEE(a, b.toExplicit2D().X(), b.toExplicit2D().Y(), c.toExplicit2D().X(), c.toExplicit2D().Y()); }
inline int orient2d_IIE(const genericPoint& a, const genericPoint& b, const genericPoint& c) { return orient2d_indirect_IIE(a, b, c.toExplicit2D().X(), c.toExplicit2D().Y()); }
inline int orient2d_III(const genericPoint& a, const genericPoint& b, const genericPoint& c) { return orient2d_indirect_III(a, b, c); }

inline int genericPoint::orient2D(const genericPoint& a, const genericPoint& b, const genericPoint& c)
{
	// Here we implicitly assume that points are 2D. Do not check.

	if (a.isExplicit2D() && b.isExplicit2D() && c.isExplicit2D()) return orient2d_EEE(a, b, c);
	if (a.isSSI() && b.isExplicit2D() && c.isExplicit2D()) return orient2d_IEE(a, b, c);
	if (a.isExplicit2D() && b.isSSI() && c.isExplicit2D()) return orient2d_IEE(b, c, a);
	if (a.isExplicit2D() && b.isExplicit2D() && c.isSSI()) return orient2d_IEE(c, a, b);
	if (a.isSSI() && b.isSSI() && c.isExplicit2D()) return orient2d_IIE(a, b, c);
	if (a.isExplicit2D() && b.isSSI() && c.isSSI())	return orient2d_IIE(b, c, a);
	if (a.isSSI() && b.isExplicit2D() && c.isSSI())	return orient2d_IIE(c, a, b);
	return orient2d_III(a, b, c);
}


inline int orient2dxy_EEE(const genericPoint& a, const genericPoint& b, const genericPoint& c) { return orient2d(a.toExplicit3D().X(), a.toExplicit3D().Y(), b.toExplicit3D().X(), b.toExplicit3D().Y(), c.toExplicit3D().X(), c.toExplicit3D().Y()); }
inline int orient2dxy_IEE(const genericPoint& a, const genericPoint& b, const genericPoint& c) { return orient2dxy_indirect_IEE(a, b.toExplicit3D().X(), b.toExplicit3D().Y(), c.toExplicit3D().X(), c.toExplicit3D().Y()); }
inline int orient2dxy_IIE(const genericPoint& a, const genericPoint& b, const genericPoint& c) { return orient2dxy_indirect_IIE(a, b, c.toExplicit3D().X(), c.toExplicit3D().Y()); }
inline int orient2dxy_III(const genericPoint& a, const genericPoint& b, const genericPoint& c) { return orient2dxy_indirect_III(a, b, c); }

inline int genericPoint::orient2Dxy(const genericPoint& a, const genericPoint& b, const genericPoint& c)
{
	if (a.isExplicit3D() && b.isExplicit3D() && c.isExplicit3D()) return orient2dxy_EEE(a, b, c);

	if (!a.isExplicit3D() && b.isExplicit3D() && c.isExplicit3D()) return orient2dxy_IEE(a, b, c);
	if (a.isExplicit3D() && !b.isExplicit3D() && c.isExplicit3D()) return orient2dxy_IEE(b, c, a);
	if (a.isExplicit3D() && b.isExplicit3D() && !c.isExplicit3D()) return orient2dxy_IEE(c, a, b);

	if (!a.isExplicit3D() && !b.isExplicit3D() && c.isExplicit3D()) return orient2dxy_IIE(a, b, c);
	if (!a.isExplicit3D() && b.isExplicit3D() && !c.isExplicit3D()) return orient2dxy_IIE(c, a, b);
	if (a.isExplicit3D() && !b.isExplicit3D() && !c.isExplicit3D()) return orient2dxy_IIE(b, c, a);

	return orient2dxy_III(a, b, c);
}


inline int orient2dyz_EEE(const genericPoint& a, const genericPoint& b, const genericPoint& c) { return orient2d(a.toExplicit3D().Y(), a.toExplicit3D().Z(), b.toExplicit3D().Y(), b.toExplicit3D().Z(), c.toExplicit3D().Y(), c.toExplicit3D().Z()); }
inline int orient2dyz_IEE(const genericPoint& a, const genericPoint& b, const genericPoint& c) { return orient2dyz_indirect_IEE(a, b.toExplicit3D().Y(), b.toExplicit3D().Z(), c.toExplicit3D().Y(), c.toExplicit3D().Z()); }
inline int orient2dyz_IIE(const genericPoint& a, const genericPoint& b, const genericPoint& c) { return orient2dyz_indirect_IIE(a, b, c.toExplicit3D().Y(), c.toExplicit3D().Z()); }
inline int orient2dyz_III(const genericPoint& a, const genericPoint& b, const genericPoint& c) { return orient2dyz_indirect_III(a, b, c); }

inline int genericPoint::orient2Dyz(const genericPoint& a, const genericPoint& b, const genericPoint& c)
{
	if (a.isExplicit3D() && b.isExplicit3D() && c.isExplicit3D()) return orient2dyz_EEE(a, b, c);

	if (!a.isExplicit3D() && b.isExplicit3D() && c.isExplicit3D()) return orient2dyz_IEE(a, b, c);
	if (a.isExplicit3D() && !b.isExplicit3D() && c.isExplicit3D()) return orient2dyz_IEE(b, c, a);
	if (a.isExplicit3D() && b.isExplicit3D() && !c.isExplicit3D()) return orient2dyz_IEE(c, a, b);

	if (!a.isExplicit3D() && !b.isExplicit3D() && c.isExplicit3D()) return orient2dyz_IIE(a, b, c);
	if (!a.isExplicit3D() && b.isExplicit3D() && !c.isExplicit3D()) return orient2dyz_IIE(c, a, b);
	if (a.isExplicit3D() && !b.isExplicit3D() && !c.isExplicit3D()) return orient2dyz_IIE(b, c, a);

	return orient2dyz_III(a, b, c);
}


inline int orient2dzx_EEE(const genericPoint& a, const genericPoint& b, const genericPoint& c) { return orient2d(a.toExplicit3D().Z(), a.toExplicit3D().X(), b.toExplicit3D().Z(), b.toExplicit3D().X(), c.toExplicit3D().Z(), c.toExplicit3D().X()); }
inline int orient2dzx_IEE(const genericPoint& a, const genericPoint& b, const genericPoint& c) { return orient2dzx_indirect_IEE(a, b.toExplicit3D().Z(), b.toExplicit3D().X(), c.toExplicit3D().Z(), c.toExplicit3D().X()); }
inline int orient2dzx_IIE(const genericPoint& a, const genericPoint& b, const genericPoint& c) { return orient2dzx_indirect_IIE(a, b, c.toExplicit3D().Z(), c.toExplicit3D().X()); }
inline int orient2dzx_III(const genericPoint& a, const genericPoint& b, const genericPoint& c) { return orient2dzx_indirect_III(a, b, c); }

inline int genericPoint::orient2Dzx(const genericPoint& a, const genericPoint& b, const genericPoint& c)
{
	if (a.isExplicit3D() && b.isExplicit3D() && c.isExplicit3D()) return orient2dzx_EEE(a, b, c);

	if (!a.isExplicit3D() && b.isExplicit3D() && c.isExplicit3D()) return orient2dzx_IEE(a, b, c);
	if (a.isExplicit3D() && !b.isExplicit3D() && c.isExplicit3D()) return orient2dzx_IEE(b, c, a);
	if (a.isExplicit3D() && b.isExplicit3D() && !c.isExplicit3D()) return orient2dzx_IEE(c, a, b);

	if (!a.isExplicit3D() && !b.isExplicit3D() && c.isExplicit3D()) return orient2dzx_IIE(a, b, c);
	if (!a.isExplicit3D() && b.isExplicit3D() && !c.isExplicit3D()) return orient2dzx_IIE(c, a, b);
	if (a.isExplicit3D() && !b.isExplicit3D() && !c.isExplicit3D()) return orient2dzx_IIE(b, c, a);

	return orient2dzx_III(a, b, c);
}


inline int orient3d_EEEE(const genericPoint& a, const genericPoint& b, const genericPoint& c, const genericPoint& d) 
{ 
	return orient3d(a.toExplicit3D().X(), a.toExplicit3D().Y(), a.toExplicit3D().Z(), b.toExplicit3D().X(), b.toExplicit3D().Y(), b.toExplicit3D().Z(),
	c.toExplicit3D().X(), c.toExplicit3D().Y(), c.toExplicit3D().Z(), d.toExplicit3D().X(), d.toExplicit3D().Y(), d.toExplicit3D().Z());
}

inline int orient3d_IEEE(const genericPoint& a, const genericPoint& b, const genericPoint& c, const genericPoint& d)
{
	return orient3d_indirect_IEEE(a, b.toExplicit3D().X(), b.toExplicit3D().Y(), b.toExplicit3D().Z(),
		c.toExplicit3D().X(), c.toExplicit3D().Y(), c.toExplicit3D().Z(), d.toExplicit3D().X(), d.toExplicit3D().Y(), d.toExplicit3D().Z());
}

inline int orient3d_IIEE(const genericPoint& a, const genericPoint& b, const genericPoint& c, const genericPoint& d)
{
	return orient3d_indirect_IIEE(a, b,
		c.toExplicit3D().X(), c.toExplicit3D().Y(), c.toExplicit3D().Z(), d.toExplicit3D().X(), d.toExplicit3D().Y(), d.toExplicit3D().Z());
}

inline int orient3d_IIIE(const genericPoint& a, const genericPoint& b, const genericPoint& c, const genericPoint& d)
{
	return orient3d_indirect_IIIE(a, b, c, d.toExplicit3D().X(), d.toExplicit3D().Y(), d.toExplicit3D().Z());
}


inline int genericPoint::orient3D(const genericPoint& a, const genericPoint& b, const genericPoint& c, const genericPoint& d)
{
	// Here we implicitly assume that points are 3D. Do not check.

	const int i = a.isExplicit3D() + b.isExplicit3D() + c.isExplicit3D() + d.isExplicit3D();

	if (i == 4) return orient3d_EEEE(a, b, c, d);
	
	if (i == 3)
	{
		if (!a.isExplicit3D()) return orient3d_IEEE(a, b, c, d);
		if (!b.isExplicit3D()) return orient3d_IEEE(b, c, a, d);
		if (!c.isExplicit3D()) return orient3d_IEEE(c, d, a, b);
		return orient3d_IEEE(d, a, c, b);
	}

	if (i == 2)
	{
		if (c.isExplicit3D() && d.isExplicit3D()) return orient3d_IIEE(a, b, c, d);
		if (b.isExplicit3D() && d.isExplicit3D()) return orient3d_IIEE(a, c, d, b);
		if (a.isExplicit3D() && d.isExplicit3D()) return orient3d_IIEE(b, c, a, d);
		if (b.isExplicit3D() && c.isExplicit3D()) return orient3d_IIEE(d, a, c, b);
		if (a.isExplicit3D() && c.isExplicit3D()) return orient3d_IIEE(d, b, a, c);
		return orient3d_IIEE(c, d, a, b);
	}

	if (i == 1)
	{
		if (d.isExplicit3D()) return orient3d_IIIE(a, b, c, d);
		if (c.isExplicit3D()) return orient3d_IIIE(d, b, a, c);
		if (b.isExplicit3D()) return orient3d_IIIE(a, c, d, b);
		return orient3d_IIIE(b, d, c, a);
	}

	return orient3d_indirect_IIII(a, b, c, d);
}

inline int inSphere_IEEEE(const genericPoint& a, const genericPoint& b, const genericPoint& c, const genericPoint& d, const genericPoint& e) {
	return inSphere_IEEEE(a,
		b.toExplicit3D().X(), b.toExplicit3D().Y(), b.toExplicit3D().Z(),
		c.toExplicit3D().X(), c.toExplicit3D().Y(), c.toExplicit3D().Z(),
		d.toExplicit3D().X(), d.toExplicit3D().Y(), d.toExplicit3D().Z(),
		e.toExplicit3D().X(), e.toExplicit3D().Y(), e.toExplicit3D().Z());
}

inline int inSphere_IIEEE(const genericPoint& a, const genericPoint& b, const genericPoint& c, const genericPoint& d, const genericPoint& e) {
	return inSphere_IIEEE(a, b,
		c.toExplicit3D().X(), c.toExplicit3D().Y(), c.toExplicit3D().Z(),
		d.toExplicit3D().X(), d.toExplicit3D().Y(), d.toExplicit3D().Z(),
		e.toExplicit3D().X(), e.toExplicit3D().Y(), e.toExplicit3D().Z());
}

inline int inSphere_IIIEE(const genericPoint& a, const genericPoint& b, const genericPoint& c, const genericPoint& d, const genericPoint& e) {
	return inSphere_IIIEE(a, b, c,
		d.toExplicit3D().X(), d.toExplicit3D().Y(), d.toExplicit3D().Z(),
		e.toExplicit3D().X(), e.toExplicit3D().Y(), e.toExplicit3D().Z());
}

inline int inSphere_IIIIE(const genericPoint& a, const genericPoint& b, const genericPoint& c, const genericPoint& d, const genericPoint& e) {
	return inSphere_IIIIE(a, b, c, d,
		e.toExplicit3D().X(), e.toExplicit3D().Y(), e.toExplicit3D().Z());
}

#define USE_LOOKUP

#ifndef USE_LOOKUP

inline int genericPoint::inSphere(const genericPoint& a, const genericPoint& b, const genericPoint& c, const genericPoint& d, const genericPoint& e)
{
	const int num_explicit = a.isExplicit3D() + b.isExplicit3D() + c.isExplicit3D() + d.isExplicit3D() + e.isExplicit3D();

	if (num_explicit == 5) return ::inSphere(a.toExplicit3D().X(), a.toExplicit3D().Y(), a.toExplicit3D().Z(), 
											b.toExplicit3D().X(), b.toExplicit3D().Y(), b.toExplicit3D().Z(),
											c.toExplicit3D().X(), c.toExplicit3D().Y(), c.toExplicit3D().Z(),
											d.toExplicit3D().X(), d.toExplicit3D().Y(), d.toExplicit3D().Z(),
											e.toExplicit3D().X(), e.toExplicit3D().Y(), e.toExplicit3D().Z());

	const genericPoint* A[5] = { &a, &b, &c, &d, &e };

	// Sort points so that I < E
	bool swapped = true;
	int sign_swap = 1;

	while (swapped) {
		swapped = false;
		for (int i = 0; i < 4; i++) {
			if (A[i]->isExplicit3D() && !A[i + 1]->isExplicit3D()) {
				std::swap(A[i], A[i + 1]);
				swapped = true;
				sign_swap *= -1;
			}
		}
	}

	if (num_explicit == 4) return sign_swap * inSphere_IEEEE(*A[0], *A[1], *A[2], *A[3], *A[4]);
	if (num_explicit == 3) return sign_swap * inSphere_IIEEE(*A[0], *A[1], *A[2], *A[3], *A[4]);
	if (num_explicit == 2) return sign_swap * inSphere_IIIEE(*A[0], *A[1], *A[2], *A[3], *A[4]);
	if (num_explicit == 1) return sign_swap * inSphere_IIIIE(*A[0], *A[1], *A[2], *A[3], *A[4]);
	return sign_swap * inSphere_IIIII(*A[0], *A[1], *A[2], *A[3], *A[4]);
}
#else

inline int genericPoint::inSphere(const genericPoint& a, const genericPoint& b, const genericPoint& c, const genericPoint& d, const genericPoint& e)
{
	const int num_explicit = a.isExplicit3D() + b.isExplicit3D() + c.isExplicit3D() + d.isExplicit3D() + e.isExplicit3D();

	if (num_explicit == 5) return ::inSphere(a.toExplicit3D().X(), a.toExplicit3D().Y(), a.toExplicit3D().Z(),
		b.toExplicit3D().X(), b.toExplicit3D().Y(), b.toExplicit3D().Z(),
		c.toExplicit3D().X(), c.toExplicit3D().Y(), c.toExplicit3D().Z(),
		d.toExplicit3D().X(), d.toExplicit3D().Y(), d.toExplicit3D().Z(),
		e.toExplicit3D().X(), e.toExplicit3D().Y(), e.toExplicit3D().Z());

	const genericPoint* A[5] = { &a, &b, &c, &d, &e };

	static const int is_lookup[] = {
		1, 0, 1, 2, 3, 4, // 00000
		1, 0, 1, 2, 3, 4, // 00001
		-1, 0, 1, 2, 4, 3, // 00010
		1, 0, 1, 2, 3, 4, // 00011
		1, 0, 1, 3, 4, 2, // 00100
		-1, 0, 1, 3, 2, 4, // 00101
		1, 0, 1, 4, 2, 3, // 00110
		1, 0, 1, 2, 3, 4, // 00111
		-1, 0, 2, 3, 4, 1, // 01000
		1, 0, 2, 3, 1, 4, // 01001
		-1, 0, 2, 4, 1, 3, // 01010
		-1, 0, 2, 1, 3, 4, // 01011
		1, 0, 3, 4, 1, 2, // 01100
		1, 0, 3, 1, 2, 4, // 01101
		-1, 0, 4, 1, 2, 3, // 01110
		1, 0, 1, 2, 3, 4, // 01111
		1, 1, 2, 3, 4, 0, // 10000
		-1, 1, 2, 3, 0, 4, // 10001
		1, 1, 2, 4, 0, 3, // 10010
		1, 1, 2, 0, 3, 4, // 10011
		-1, 1, 3, 4, 0, 2, // 10100
		-1, 1, 3, 0, 2, 4, // 10101
		1, 1, 4, 0, 2, 3, // 10110
		-1, 1, 0, 2, 3, 4, // 10111
		1, 2, 3, 4, 0, 1, // 11000
		1, 2, 3, 0, 1, 4, // 11001
		-1, 2, 4, 0, 1, 3, // 11010
		1, 2, 0, 1, 3, 4, // 11011
		1, 3, 4, 0, 1, 2, // 11100
		-1, 3, 0, 1, 2, 4, // 11101
		1, 4, 0, 1, 2, 3, // 11110
		1, 0, 1, 2, 3, 4 // 11111
	};

	const int idx = (a.isExplicit3D() << 4) + (b.isExplicit3D() << 3) + (c.isExplicit3D() << 2) + (d.isExplicit3D() << 1) + e.isExplicit3D();
	const int* cfg = is_lookup + idx * 6;

	if (num_explicit == 4) return cfg[0] * inSphere_IEEEE(*A[cfg[1]], *A[cfg[2]], *A[cfg[3]], *A[cfg[4]], *A[cfg[5]]);
	if (num_explicit == 3) return cfg[0] * inSphere_IIEEE(*A[cfg[1]], *A[cfg[2]], *A[cfg[3]], *A[cfg[4]], *A[cfg[5]]);
	if (num_explicit == 2) return cfg[0] * inSphere_IIIEE(*A[cfg[1]], *A[cfg[2]], *A[cfg[3]], *A[cfg[4]], *A[cfg[5]]);
	if (num_explicit == 1) return cfg[0] * inSphere_IIIIE(*A[cfg[1]], *A[cfg[2]], *A[cfg[3]], *A[cfg[4]], *A[cfg[5]]);
	return cfg[0] * inSphere_IIIII(*A[cfg[1]], *A[cfg[2]], *A[cfg[3]], *A[cfg[4]], *A[cfg[5]]);
}
#endif


inline int inGabrielSphere_EEEE(const genericPoint& a, const genericPoint& b, const genericPoint& c, const genericPoint& d)
{
	return inGabrielSphere(a.toExplicit3D().X(), a.toExplicit3D().Y(), a.toExplicit3D().Z(), b.toExplicit3D().X(), b.toExplicit3D().Y(), b.toExplicit3D().Z(),
		c.toExplicit3D().X(), c.toExplicit3D().Y(), c.toExplicit3D().Z(), d.toExplicit3D().X(), d.toExplicit3D().Y(), d.toExplicit3D().Z());
}

inline int inGabrielSphere_IEEE(const genericPoint& a, const genericPoint& b, const genericPoint& c, const genericPoint& d)
{
	return inGabrielSphere_IEEE(a, b.toExplicit3D().X(), b.toExplicit3D().Y(), b.toExplicit3D().Z(),
		c.toExplicit3D().X(), c.toExplicit3D().Y(), c.toExplicit3D().Z(), d.toExplicit3D().X(), d.toExplicit3D().Y(), d.toExplicit3D().Z());
}

inline int inGabrielSphere_EIEE(const genericPoint& a, const genericPoint& b, const genericPoint& c, const genericPoint& d)
{
	return inGabrielSphere_EIEE(b, a.toExplicit3D().X(), a.toExplicit3D().Y(), a.toExplicit3D().Z(),
		c.toExplicit3D().X(), c.toExplicit3D().Y(), c.toExplicit3D().Z(), d.toExplicit3D().X(), d.toExplicit3D().Y(), d.toExplicit3D().Z());
}

inline int inGabrielSphere_IIEE(const genericPoint& a, const genericPoint& b, const genericPoint& c, const genericPoint& d)
{
	return inGabrielSphere_IIEE(a, b,
		c.toExplicit3D().X(), c.toExplicit3D().Y(), c.toExplicit3D().Z(), d.toExplicit3D().X(), d.toExplicit3D().Y(), d.toExplicit3D().Z());
}

inline int inGabrielSphere_EIIE(const genericPoint& a, const genericPoint& b, const genericPoint& c, const genericPoint& d)
{
	return inGabrielSphere_EIIE(b, c,
		a.toExplicit3D().X(), a.toExplicit3D().Y(), a.toExplicit3D().Z(), d.toExplicit3D().X(), d.toExplicit3D().Y(), d.toExplicit3D().Z());
}

inline int inGabrielSphere_IIIE(const genericPoint& a, const genericPoint& b, const genericPoint& c, const genericPoint& d)
{
	return inGabrielSphere_IIIE(a, b, c, d.toExplicit3D().X(), d.toExplicit3D().Y(), d.toExplicit3D().Z());
}

inline int inGabrielSphere_EIII(const genericPoint& a, const genericPoint& b, const genericPoint& c, const genericPoint& d)
{
	return inGabrielSphere_EIII(b, c, d, a.toExplicit3D().X(), a.toExplicit3D().Y(), a.toExplicit3D().Z());
}


inline int genericPoint::inGabrielSphere(const genericPoint& a, const genericPoint& b, const genericPoint& c, const genericPoint& d)
{
	const int i = a.isExplicit3D() + b.isExplicit3D() + c.isExplicit3D() + d.isExplicit3D();

	if (i == 4) return inGabrielSphere_EEEE(a, b, c, d);

	if (i == 3)
	{
		if (!a.isExplicit3D()) return inGabrielSphere_IEEE(a, b, c, d);
		if (!b.isExplicit3D()) return inGabrielSphere_EIEE(a, b, c, d);
		if (!c.isExplicit3D()) return inGabrielSphere_EIEE(a, c, b, d);
		return inGabrielSphere_EIEE(a, d, b, c);
	}

	if (i == 2)
	{
		if (c.isExplicit3D() && d.isExplicit3D()) return inGabrielSphere_IIEE(a, b, c, d);
		if (b.isExplicit3D() && d.isExplicit3D()) return inGabrielSphere_IIEE(a, c, b, d);
		if (a.isExplicit3D() && d.isExplicit3D()) return inGabrielSphere_EIIE(a, b, c, d);
		if (b.isExplicit3D() && c.isExplicit3D()) return inGabrielSphere_IIEE(a, d, b, c);
		if (a.isExplicit3D() && c.isExplicit3D()) return inGabrielSphere_EIIE(a, b, d, c);
		return inGabrielSphere_EIIE(a, c, d, b);
	}

	if (i == 1)
	{
		if (d.isExplicit3D()) return inGabrielSphere_IIIE(a, b, c, d);
		if (c.isExplicit3D()) return inGabrielSphere_IIIE(a, b, d, c);
		if (b.isExplicit3D()) return inGabrielSphere_IIIE(a, c, d, b);
		return inGabrielSphere_EIII(a, b, c, d);
	}

	return inGabrielSphere_IIII(a, b, c, d);
}


inline int incircle2d_EEEE(const genericPoint& a, const genericPoint& b, const genericPoint& c, const genericPoint& d)
{
	return incircle(a.toExplicit2D().X(), a.toExplicit2D().Y(), b.toExplicit2D().X(), b.toExplicit2D().Y(), c.toExplicit2D().X(), c.toExplicit2D().Y(), d.toExplicit2D().X(), d.toExplicit2D().Y());
}

inline int incircle2d_IEEE(const genericPoint& a, const genericPoint& b, const genericPoint& c, const genericPoint& d)
{
	return incircle_indirect_IEEE(a, b.toExplicit2D().X(), b.toExplicit2D().Y(), c.toExplicit2D().X(), c.toExplicit2D().Y(), d.toExplicit2D().X(), d.toExplicit2D().Y());
}

inline int incircle2d_IIEE(const genericPoint& a, const genericPoint& b, const genericPoint& c, const genericPoint& d)
{
	return incircle_indirect_IIEE(a, b, c.toExplicit2D().X(), c.toExplicit2D().Y(), d.toExplicit2D().X(), d.toExplicit2D().Y());
}

inline int incircle2d_IIIE(const genericPoint& a, const genericPoint& b, const genericPoint& c, const genericPoint& d)
{
	return incircle_indirect_IIIE(a, b, c, d.toExplicit2D().X(), d.toExplicit2D().Y());
}


inline int genericPoint::incircle(const genericPoint& a, const genericPoint& b, const genericPoint& c, const genericPoint& d)
{
	const int i = a.isExplicit2D() + b.isExplicit2D() + c.isExplicit2D() + d.isExplicit2D();

	if (i == 4) return incircle2d_EEEE(a, b, c, d);

	if (i == 3)
	{
		if (!a.isExplicit2D()) return incircle2d_IEEE(a, b, c, d);
		if (!b.isExplicit2D()) return incircle2d_IEEE(b, c, a, d);
		if (!c.isExplicit2D()) return incircle2d_IEEE(c, d, a, b);
		return incircle2d_IEEE(d, a, c, b);
	}

	if (i == 2)
	{
		if (c.isExplicit2D() && d.isExplicit2D()) return incircle2d_IIEE(a, b, c, d);
		if (b.isExplicit2D() && d.isExplicit2D()) return incircle2d_IIEE(a, c, d, b);
		if (a.isExplicit2D() && d.isExplicit2D()) return incircle2d_IIEE(b, c, a, d);
		if (b.isExplicit2D() && c.isExplicit2D()) return incircle2d_IIEE(d, a, c, b);
		if (a.isExplicit2D() && c.isExplicit2D()) return incircle2d_IIEE(d, b, a, c);
		return incircle2d_IIEE(c, d, a, b);
	}

	if (i == 1)
	{
		if (d.isExplicit2D()) return incircle2d_IIIE(a, b, c, d);
		if (c.isExplicit2D()) return incircle2d_IIIE(d, b, a, c);
		if (b.isExplicit2D()) return incircle2d_IIIE(a, c, d, b);
		return incircle2d_IIIE(b, d, c, a);
	}

	return incircle_indirect_IIII(a, b, c, d);
}

inline int incircle2dxy_EEEE(const genericPoint& a, const genericPoint& b, const genericPoint& c, const genericPoint& d)
{
	return incircle(a.toExplicit3D().X(), a.toExplicit3D().Y(), b.toExplicit3D().X(), b.toExplicit3D().Y(), c.toExplicit3D().X(), c.toExplicit3D().Y(), d.toExplicit3D().X(), d.toExplicit3D().Y());
}

inline int incircle2dxy_IEEE(const genericPoint& a, const genericPoint& b, const genericPoint& c, const genericPoint& d)
{
	return incirclexy_indirect_IEEE(a, b.toExplicit3D().X(), b.toExplicit3D().Y(), c.toExplicit3D().X(), c.toExplicit3D().Y(), d.toExplicit3D().X(), d.toExplicit3D().Y());
}

inline int incircle2dxy_IIEE(const genericPoint& a, const genericPoint& b, const genericPoint& c, const genericPoint& d)
{
	return incirclexy_indirect_IIEE(a, b, c.toExplicit3D().X(), c.toExplicit3D().Y(), d.toExplicit3D().X(), d.toExplicit3D().Y());
}

inline int incircle2dxy_IIIE(const genericPoint& a, const genericPoint& b, const genericPoint& c, const genericPoint& d)
{
	return incirclexy_indirect_IIIE(a, b, c, d.toExplicit3D().X(), d.toExplicit3D().Y());
}


inline int genericPoint::incirclexy(const genericPoint& a, const genericPoint& b, const genericPoint& c, const genericPoint& d)
{
	int i = a.isExplicit3D() + b.isExplicit3D() + c.isExplicit3D() + d.isExplicit3D();

	if (i == 4) return incircle2dxy_EEEE(a, b, c, d);

	if (i == 3)
	{
		if (!a.isExplicit3D()) return incircle2dxy_IEEE(a, b, c, d);
		if (!b.isExplicit3D()) return incircle2dxy_IEEE(b, c, a, d);
		if (!c.isExplicit3D()) return incircle2dxy_IEEE(c, d, a, b);
		return incircle2dxy_IEEE(d, a, c, b);
	}

	if (i == 2)
	{
		if (c.isExplicit3D() && d.isExplicit3D()) return incircle2dxy_IIEE(a, b, c, d);
		if (b.isExplicit3D() && d.isExplicit3D()) return incircle2dxy_IIEE(a, c, d, b);
		if (a.isExplicit3D() && d.isExplicit3D()) return incircle2dxy_IIEE(b, c, a, d);
		if (b.isExplicit3D() && c.isExplicit3D()) return incircle2dxy_IIEE(d, a, c, b);
		if (a.isExplicit3D() && c.isExplicit3D()) return incircle2dxy_IIEE(d, b, a, c);
		return incircle2dxy_IIEE(c, d, a, b);
	}

	if (i == 1)
	{
		if (d.isExplicit3D()) return incircle2dxy_IIIE(a, b, c, d);
		if (c.isExplicit3D()) return incircle2dxy_IIIE(d, b, a, c);
		if (b.isExplicit3D()) return incircle2dxy_IIIE(a, c, d, b);
		return incircle2dxy_IIIE(b, d, c, a);
	}

	return incirclexy_indirect_IIII(a, b, c, d);
}


// These functions assume that point is an SSI
inline bool genericPoint::getIntervalLambda(interval_number& lx, interval_number& ly, interval_number& d) const {
	if (isExplicit2D()) {
		const explicitPoint2D& e = toExplicit2D();
		lx = e.X(), ly = e.Y(), d = 1;
		return true;
	}
	else {
		assert(isSSI());
		return toSSI().getIntervalLambda(lx, ly, d);
	}
}

inline void genericPoint::getBigfloatLambda(bigfloat& lx, bigfloat& ly, bigfloat& d) const {
	if (isExplicit2D()) {
		const explicitPoint2D& e = toExplicit2D();
		lx = e.X(), ly = e.Y(), d = 1;
	}
	else {
		assert(isSSI());
		toSSI().getBigfloatLambda(lx, ly, d);
	}
}

// These functions assume that point is an implicit 3D
inline bool genericPoint::getIntervalLambda(interval_number& lx, interval_number& ly, interval_number& lz, interval_number& d) const {
	if (isExplicit3D()) {
		const explicitPoint3D& e = toExplicit3D();
		lx = e.X(), ly = e.Y(), lz = e.Z(); d = 1;
		return true;
	}
	else if (isLPI()) return toLPI().getIntervalLambda(lx, ly, lz, d);
	else if (isTPI()) return toTPI().getIntervalLambda(lx, ly, lz, d);
	else if (isLNC()) return toLNC().getIntervalLambda(lx, ly, lz, d);
	else if (isBPT()) return toBPT().getIntervalLambda(lx, ly, lz, d);
	else {
		assert(isTBC());
		return toTBC().getIntervalLambda(lx, ly, lz, d);
	}
}

inline void genericPoint::getBigfloatLambda(bigfloat& lx, bigfloat& ly, bigfloat& lz, bigfloat& d) const {
	if (isExplicit3D()) {
		const explicitPoint3D& e = toExplicit3D();
		lx = e.X(), ly = e.Y(), lz = e.Z(); d = 1;
	}
	else if (isLPI()) toLPI().getBigfloatLambda(lx, ly, lz, d);
	else if (isTPI()) toTPI().getBigfloatLambda(lx, ly, lz, d);
	else if (isLNC()) toLNC().getBigfloatLambda(lx, ly, lz, d);
	else if (isBPT()) toBPT().getBigfloatLambda(lx, ly, lz, d);
	else {
		assert(isTBC());
		toTBC().getBigfloatLambda(lx, ly, lz, d);
	}
}

// Type-specific lambdas

inline bool implicitPoint2D_SSI::getIntervalLambda(interval_number& lx, interval_number& ly, interval_number &d) const
{
	lx = dfilter_lambda_x;
	ly = dfilter_lambda_y;
	d = dfilter_denominator;
	return (dfilter_denominator.signIsReliable());
}

inline implicitPoint2D_SSI::implicitPoint2D_SSI(const genericPoint& l11, const genericPoint& l12,
	const genericPoint& l21, const genericPoint& l22)
	: genericPoint(Point_Type::SSI), l1_1(l11), l1_2(l12), l2_1(l21), l2_2(l22)
{
	interval_number l1[3], l2[3], m1[3], m2[3];
	l1_1.getIntervalLambda(l1[0], l1[1], l1[2]);
	l1_2.getIntervalLambda(l2[0], l2[1], l2[2]);
	l2_1.getIntervalLambda(m1[0], m1[1], m1[2]);
	l2_2.getIntervalLambda(m2[0], m2[1], m2[2]);
	lambda2d_SSI<interval_number, interval_number>(l1, l2, m1, m2, dfilter_lambda_x, dfilter_lambda_y, dfilter_denominator);
	if (dfilter_denominator.isNegative()) {
		dfilter_lambda_x.negate();
		dfilter_lambda_y.negate();
		dfilter_denominator.negate();
	}
}

inline bool implicitPoint3D_LPI::getIntervalLambda(interval_number& lx, interval_number& ly, interval_number& lz, interval_number &d) const
{
	lx = dfilter_lambda_x;
	ly = dfilter_lambda_y;
	lz = dfilter_lambda_z;
	d = dfilter_denominator;
	return (dfilter_denominator.signIsReliable());
}

inline implicitPoint3D_LPI::implicitPoint3D_LPI(const genericPoint& _p, const genericPoint& _q,
	const genericPoint& _r, const genericPoint& _s, const genericPoint& _t)
	: genericPoint(Point_Type::LPI), ip(_p), iq(_q), ir(_r), is(_s), it(_t)
{
	interval_number p[4], q[4], r[4], s[4], t[4];
	ip.getIntervalLambda(p[0], p[1], p[2], p[3]);
	iq.getIntervalLambda(q[0], q[1], q[2], q[3]);
	ir.getIntervalLambda(r[0], r[1], r[2], r[3]);
	is.getIntervalLambda(s[0], s[1], s[2], s[3]);
	it.getIntervalLambda(t[0], t[1], t[2], t[3]);

	lambda3d_LPI<interval_number, interval_number>(p, q, r, s, t, dfilter_lambda_x, dfilter_lambda_y, dfilter_lambda_z, dfilter_denominator);
	if (dfilter_denominator.isNegative()) {
		dfilter_lambda_x.negate();
		dfilter_lambda_y.negate();
		dfilter_lambda_z.negate();
		dfilter_denominator.negate();
	}
}

inline bool implicitPoint3D_TPI::getIntervalLambda(interval_number& lx, interval_number& ly, interval_number& lz, interval_number& d) const
{
	lx = dfilter_lambda_x;
	ly = dfilter_lambda_y;
	lz = dfilter_lambda_z;
	d = dfilter_denominator;
	return (dfilter_denominator.signIsReliable());
}

inline implicitPoint3D_TPI::implicitPoint3D_TPI(const genericPoint& _v1, const genericPoint& _v2, const genericPoint& _v3,
	const genericPoint& _w1, const genericPoint& _w2, const genericPoint& _w3,
	const genericPoint& _u1, const genericPoint& _u2, const genericPoint& _u3)
	: genericPoint(Point_Type::TPI), iv1(_v1), iv2(_v2), iv3(_v3), iw1(_w1), iw2(_w2), iw3(_w3), iu1(_u1), iu2(_u2), iu3(_u3)
{
	interval_number v1[4], v2[4], v3[4], w1[4], w2[4], w3[4], u1[4], u2[4], u3[4];
	iv1.getIntervalLambda(v1[0], v1[1], v1[2], v1[3]);
	iv2.getIntervalLambda(v2[0], v2[1], v2[2], v2[3]);
	iv3.getIntervalLambda(v3[0], v3[1], v3[2], v3[3]);
	iw1.getIntervalLambda(w1[0], w1[1], w1[2], w1[3]);
	iw2.getIntervalLambda(w2[0], w2[1], w2[2], w2[3]);
	iw3.getIntervalLambda(w3[0], w3[1], w3[2], w3[3]);
	iu1.getIntervalLambda(u1[0], u1[1], u1[2], u1[3]);
	iu2.getIntervalLambda(u2[0], u2[1], u2[2], u2[3]);
	iu3.getIntervalLambda(u3[0], u3[1], u3[2], u3[3]);

	lambda3d_TPI<interval_number, interval_number>(
		v1, v2, v3, w1, w2, w3, u1 ,u2, u3,
		dfilter_lambda_x, dfilter_lambda_y, dfilter_lambda_z, dfilter_denominator);
	if (dfilter_denominator.isNegative()) {
		dfilter_lambda_x.negate();
		dfilter_lambda_y.negate();
		dfilter_lambda_z.negate();
		dfilter_denominator.negate();
	}
}

inline bool implicitPoint3D_LNC::getIntervalLambda(interval_number& lx, interval_number& ly, interval_number& lz, interval_number& d) const
{
	lx = dfilter_lambda_x;
	ly = dfilter_lambda_y;
	lz = dfilter_lambda_z;
	d = dfilter_denominator;
	return true;
}

inline implicitPoint3D_LNC::implicitPoint3D_LNC(const genericPoint& _p, const genericPoint& _q,
	const double _t)
	: genericPoint(Point_Type::LNC), ip(_p), iq(_q), t(_t)
{
	interval_number p[4], q[4];
	ip.getIntervalLambda(p[0], p[1], p[2], p[3]);
	iq.getIntervalLambda(q[0], q[1], q[2], q[3]);
	lambda3d_LNC<interval_number, interval_number, interval_number>(p, q, t, dfilter_lambda_x, dfilter_lambda_y, dfilter_lambda_z, dfilter_denominator);
}

inline bool implicitPoint3D_BPT::getIntervalLambda(interval_number& lx, interval_number& ly, interval_number& lz, interval_number& d) const
{
	lx = dfilter_lambda_x;
	ly = dfilter_lambda_y;
	lz = dfilter_lambda_z;
	d = dfilter_denominator;
	return true;
}

inline implicitPoint3D_BPT::implicitPoint3D_BPT(const genericPoint& _p, const genericPoint& _q, const genericPoint& _r,
	const double _v, const double _u)
	: genericPoint(Point_Type::BPT), ip(_p), iq(_q), ir(_r), v(_v), u(_u)
{
	interval_number p[4], q[4], r[4];
	ip.getIntervalLambda(p[0], p[1], p[2], p[3]);
	iq.getIntervalLambda(q[0], q[1], q[2], q[3]);
	ir.getIntervalLambda(r[0], r[1], r[2], r[3]);
	lambda3d_BPT<interval_number, interval_number, interval_number>(p, q, r, u, v, dfilter_lambda_x, dfilter_lambda_y, dfilter_lambda_z, dfilter_denominator);
}

inline bool implicitPoint3D_TBC::getIntervalLambda(interval_number& lx, interval_number& ly, interval_number& lz, interval_number& d) const
{
	lx = dfilter_lambda_x;
	ly = dfilter_lambda_y;
	lz = dfilter_lambda_z;
	d = dfilter_denominator;
	return true;
}

inline implicitPoint3D_TBC::implicitPoint3D_TBC(const genericPoint& _p, const genericPoint& _q, const genericPoint& _r, const genericPoint& _s)
	: genericPoint(Point_Type::TBC), ip(_p), iq(_q), ir(_r), is(_s)
{
	interval_number p[4], q[4], r[4], s[4];
	ip.getIntervalLambda(p[0], p[1], p[2], p[3]);
	iq.getIntervalLambda(q[0], q[1], q[2], q[3]);
	ir.getIntervalLambda(r[0], r[1], r[2], r[3]);
	is.getIntervalLambda(s[0], s[1], s[2], s[3]);
	lambda3d_TBC<interval_number, interval_number>(p, q, r, s, dfilter_lambda_x, dfilter_lambda_y, dfilter_lambda_z, dfilter_denominator);
}


inline void genericPoint::getExpansionLambda(expansion& lx, expansion& ly, expansion& d) const {
	if (isExplicit2D()) {
		const explicitPoint2D& e = toExplicit2D();
		lx = e.X(), ly = e.Y(), d = 1;
	}
	else {
		assert(isSSI());
		toSSI().getExpansionLambda(lx, ly, d);
	}
}

inline void genericPoint::getExpansionLambda(expansion& lx, expansion& ly, expansion& lz, expansion& d) const {
	if (isExplicit3D()) {
		const explicitPoint3D& e = toExplicit3D();
		lx = e.X(), ly = e.Y(), lz = e.Z(), d = 1;
	}
	else if (isLPI()) toLPI().getExpansionLambda(lx, ly, lz, d);
	else if (isTPI()) toTPI().getExpansionLambda(lx, ly, lz, d);
	else if (isLNC()) toLNC().getExpansionLambda(lx, ly, lz, d);
	else if (isBPT()) toBPT().getExpansionLambda(lx, ly, lz, d);
	else {
		assert(isTBC());
		toTBC().getExpansionLambda(lx, ly, lz, d);
	}
}

inline void implicitPoint2D_SSI::getExpansionLambda(expansion& lx, expansion& ly, expansion& d) const
{
	if (l1_1.isExplicit2D() && l1_2.isExplicit2D() && l2_1.isExplicit2D() && l2_2.isExplicit2D()) {
		const explicitPoint2D& p1 = l1_1.toExplicit2D(), & p2 = l1_2.toExplicit2D();
		const explicitPoint2D& q1 = l2_1.toExplicit2D(), & q2 = l2_2.toExplicit2D();
		const s_expansion l1[3] = { p1.X(), p1.Y(), 1 }, l2[3] = { p2.X(), p2.Y(), 1 };
		const s_expansion m1[3] = { q1.X(), q1.Y(), 1 }, m2[3] = { q2.X(), q2.Y(), 1 };
		lambda2d_SSI<s_expansion, expansion>(l1, l2, m1, m2, lx, ly, d);
	}
	else {
		expansion l1[3], l2[3], m1[3], m2[3];
		l1_1.getExpansionLambda(l1[0], l1[1], l1[2]);
		l1_2.getExpansionLambda(l2[0], l2[1], l2[2]);
		l2_1.getExpansionLambda(m1[0], m1[1], m1[2]);
		l2_2.getExpansionLambda(m2[0], m2[1], m2[2]);
		lambda2d_SSI<expansion, expansion>(l1, l2, m1, m2, lx, ly, d);
	}
	if (sgn(d) < 0)
	{
		lx.negate();
		ly.negate();
		d.negate();
	}
}









inline void implicitPoint3D_LPI::getExpansionLambda(expansion& lx, expansion& ly, expansion& lz, expansion& d) const
{
	if (ip.isExplicit3D() && iq.isExplicit3D() && ir.isExplicit3D() && is.isExplicit3D() && it.isExplicit3D()) {
		const explicitPoint3D& ep = ip.toExplicit3D(), & eq = iq.toExplicit3D(), & er = ir.toExplicit3D();
		const explicitPoint3D& es = is.toExplicit3D(), & et = it.toExplicit3D();
		const s_expansion p[4] = { ep.X(), ep.Y(), ep.Z(), 1 };
		const s_expansion q[4] = { eq.X(), eq.Y(), eq.Z(), 1 };
		const s_expansion r[4] = { er.X(), er.Y(), er.Z(), 1 };
		const s_expansion s[4] = { es.X(), es.Y(), es.Z(), 1 };
		const s_expansion t[4] = { et.X(), et.Y(), et.Z(), 1 };
		lambda3d_LPI<s_expansion, expansion>(p, q, r, s, t, lx, ly, lz, d);
	}
	else {
		expansion p[4], q[4], r[4], s[4], t[4];
		ip.getExpansionLambda(p[0], p[1], p[2], p[3]);
		iq.getExpansionLambda(q[0], q[1], q[2], q[3]);
		ir.getExpansionLambda(r[0], r[1], r[2], r[3]);
		is.getExpansionLambda(s[0], s[1], s[2], s[3]);
		it.getExpansionLambda(t[0], t[1], t[2], t[3]);
		lambda3d_LPI<expansion, expansion>(p, q, r, s, t, lx, ly, lz, d);
	}
	if (sgn(d) < 0)
	{
		lx.negate();
		ly.negate();
		lz.negate();
		d.negate();
	}
}

inline void implicitPoint3D_TPI::getExpansionLambda(expansion& lx, expansion& ly, expansion& lz, expansion& d) const
{
	if (iv1.isExplicit3D() && iv2.isExplicit3D() && iv3.isExplicit3D() && 
		iw1.isExplicit3D() && iw2.isExplicit3D() && iw3.isExplicit3D() &&
		iu1.isExplicit3D() && iu2.isExplicit3D() && iu3.isExplicit3D()) {
		const explicitPoint3D& ev1 = iv1.toExplicit3D(), & ev2 = iv2.toExplicit3D(), & ev3 = iv3.toExplicit3D();
		const explicitPoint3D& ew1 = iw1.toExplicit3D(), & ew2 = iw2.toExplicit3D(), & ew3 = iw3.toExplicit3D();
		const explicitPoint3D& eu1 = iu1.toExplicit3D(), & eu2 = iu2.toExplicit3D(), & eu3 = iu3.toExplicit3D();
		const s_expansion v1[4] = { ev1.X(), ev1.Y(), ev1.Z(), 1 };
		const s_expansion v2[4] = { ev2.X(), ev2.Y(), ev2.Z(), 1 };
		const s_expansion v3[4] = { ev3.X(), ev3.Y(), ev3.Z(), 1 };
		const s_expansion w1[4] = { ew1.X(), ew1.Y(), ew1.Z(), 1 };
		const s_expansion w2[4] = { ew2.X(), ew2.Y(), ew2.Z(), 1 };
		const s_expansion w3[4] = { ew3.X(), ew3.Y(), ew3.Z(), 1 };
		const s_expansion u1[4] = { eu1.X(), eu1.Y(), eu1.Z(), 1 };
		const s_expansion u2[4] = { eu2.X(), eu2.Y(), eu2.Z(), 1 };
		const s_expansion u3[4] = { eu3.X(), eu3.Y(), eu3.Z(), 1 };
		lambda3d_TPI<s_expansion, expansion>(v1, v2, v3, w1, w2, w3, u1, u2, u3, lx, ly, lz, d);
	}
	else {
		expansion v1[4], v2[4], v3[4], w1[4], w2[4], w3[4], u1[4], u2[4], u3[4];
		iv1.getExpansionLambda(v1[0], v1[1], v1[2], v1[3]);
		iv2.getExpansionLambda(v2[0], v2[1], v2[2], v2[3]);
		iv3.getExpansionLambda(v3[0], v3[1], v3[2], v3[3]);
		iw1.getExpansionLambda(w1[0], w1[1], w1[2], w1[3]);
		iw2.getExpansionLambda(w2[0], w2[1], w2[2], w2[3]);
		iw3.getExpansionLambda(w3[0], w3[1], w3[2], w3[3]);
		iu1.getExpansionLambda(u1[0], u1[1], u1[2], u1[3]);
		iu2.getExpansionLambda(u2[0], u2[1], u2[2], u2[3]);
		iu3.getExpansionLambda(u3[0], u3[1], u3[2], u3[3]);
		lambda3d_TPI<expansion, expansion>(v1, v2, v3, w1, w2, w3, u1, u2, u3, lx, ly, lz, d);
	}
	if (sgn(d) < 0)
	{
		lx.negate();
		ly.negate();
		lz.negate();
		d.negate();
	}
}

inline void implicitPoint3D_LNC::getExpansionLambda(expansion& lx, expansion& ly, expansion& lz, expansion& d) const
{
	if (ip.isExplicit3D() && iq.isExplicit3D()) {
		const explicitPoint3D& ep = ip.toExplicit3D(), & eq = iq.toExplicit3D();
		const s_expansion p[4] = { ep.X(), ep.Y(), ep.Z(), 1 };
		const s_expansion q[4] = { eq.X(), eq.Y(), eq.Z(), 1 };
		lambda3d_LNC<s_expansion, s_expansion, expansion>(p, q, T(), lx, ly, lz, d);
	}
	else {
		expansion p[4], q[4];
		ip.getExpansionLambda(p[0], p[1], p[2], p[3]);
		iq.getExpansionLambda(q[0], q[1], q[2], q[3]);
		lambda3d_LNC<expansion, s_expansion, expansion>(p, q, T(), lx, ly, lz, d);
	}
}

inline void implicitPoint3D_BPT::getExpansionLambda(expansion& lx, expansion& ly, expansion& lz, expansion& d) const
{
	if (ip.isExplicit3D() && iq.isExplicit3D() && ir.isExplicit3D()) {
		const explicitPoint3D& ep = ip.toExplicit3D(), & eq = iq.toExplicit3D(), & er = ir.toExplicit3D();
		const s_expansion p[4] = { ep.X(), ep.Y(), ep.Z(), 1 };
		const s_expansion q[4] = { eq.X(), eq.Y(), eq.Z(), 1 };
		const s_expansion r[4] = { er.X(), er.Y(), er.Z(), 1 };
		lambda3d_BPT<s_expansion, s_expansion, expansion>(p, q, r, U(), V(), lx, ly, lz, d);
	}
	else {
		expansion p[4], q[4], r[4];
		ip.getExpansionLambda(p[0], p[1], p[2], p[3]);
		iq.getExpansionLambda(q[0], q[1], q[2], q[3]);
		ir.getExpansionLambda(r[0], r[1], r[2], r[3]);
		lambda3d_BPT<expansion, s_expansion, expansion>(p, q, r, U(), V(), lx, ly, lz, d);
	}
}

inline void implicitPoint3D_TBC::getExpansionLambda(expansion& lx, expansion& ly, expansion& lz, expansion& d) const
{
	if (ip.isExplicit3D() && iq.isExplicit3D() && ir.isExplicit3D() && is.isExplicit3D()) {
		const explicitPoint3D& ep = ip.toExplicit3D(), & eq = iq.toExplicit3D(), & er = ir.toExplicit3D(), & es = is.toExplicit3D();
		const s_expansion p[4] = { ep.X(), ep.Y(), ep.Z(), 1 };
		const s_expansion q[4] = { eq.X(), eq.Y(), eq.Z(), 1 };
		const s_expansion r[4] = { er.X(), er.Y(), er.Z(), 1 };
		const s_expansion s[4] = { es.X(), es.Y(), es.Z(), 1 };
		lambda3d_TBC<s_expansion, expansion>(p, q, r, s, lx, ly, lz, d);
	}
	else {
		expansion p[4], q[4], r[4], s[4];
		ip.getExpansionLambda(p[0], p[1], p[2], p[3]);
		iq.getExpansionLambda(q[0], q[1], q[2], q[3]);
		ir.getExpansionLambda(r[0], r[1], r[2], r[3]);
		is.getExpansionLambda(s[0], s[1], s[2], s[3]);
		lambda3d_TBC<expansion, expansion>(p, q, r, s, lx, ly, lz, d);
	}
}














inline void implicitPoint2D_SSI::getBigfloatLambda(bigfloat& lx, bigfloat& ly, bigfloat& d) const
{
	bigfloat l1[3], l2[3], m1[3], m2[3];
	l1_1.getBigfloatLambda(l1[0], l1[1], l1[2]);
	l1_2.getBigfloatLambda(l2[0], l2[1], l2[2]);
	l2_1.getBigfloatLambda(m1[0], m1[1], m1[2]);
	l2_2.getBigfloatLambda(m2[0], m2[1], m2[2]);
	lambda2d_SSI<bigfloat, bigfloat>(l1, l2, m1, m2, lx, ly, d);
	if (sgn(d) < 0)
	{
		lx = -lx;
		ly = -ly;
		d = -d;
	}
}

inline void implicitPoint3D_LPI::getBigfloatLambda(bigfloat& lx, bigfloat& ly, bigfloat& lz, bigfloat& d) const
{
	bigfloat p[4], q[4], r[4], s[4], t[4];
	ip.getBigfloatLambda(p[0], p[1], p[2], p[3]);
	iq.getBigfloatLambda(q[0], q[1], q[2], q[3]);
	ir.getBigfloatLambda(r[0], r[1], r[2], r[3]);
	is.getBigfloatLambda(s[0], s[1], s[2], s[3]);
	it.getBigfloatLambda(t[0], t[1], t[2], t[3]);
	lambda3d_LPI<bigfloat, bigfloat>(p, q, r, s, t, lx, ly, lz, d);
	if (sgn(d) < 0)
	{
		lx = -lx;
		ly = -ly;
		lz = -lz;
		d = -d;
	}
}

inline void implicitPoint3D_TPI::getBigfloatLambda(bigfloat& lx, bigfloat& ly, bigfloat& lz, bigfloat& d) const
{
	bigfloat v1[4], v2[4], v3[4], w1[4], w2[4], w3[4], u1[4], u2[4], u3[4];
	iv1.getBigfloatLambda(v1[0], v1[1], v1[2], v1[3]);
	iv2.getBigfloatLambda(v2[0], v2[1], v2[2], v2[3]);
	iv3.getBigfloatLambda(v3[0], v3[1], v3[2], v3[3]);
	iw1.getBigfloatLambda(w1[0], w1[1], w1[2], w1[3]);
	iw2.getBigfloatLambda(w2[0], w2[1], w2[2], w2[3]);
	iw3.getBigfloatLambda(w3[0], w3[1], w3[2], w3[3]);
	iu1.getBigfloatLambda(u1[0], u1[1], u1[2], u1[3]);
	iu2.getBigfloatLambda(u2[0], u2[1], u2[2], u2[3]);
	iu3.getBigfloatLambda(u3[0], u3[1], u3[2], u3[3]);

	lambda3d_TPI<bigfloat, bigfloat>(v1, v2, v3, w1, w2, w3, u1, u2, u3, lx, ly, lz, d);
	if (sgn(d) < 0)
	{
		lx = -lx;
		ly = -ly;
		lz = -lz;
		d = -d;
	}
}

inline void implicitPoint3D_LNC::getBigfloatLambda(bigfloat& lx, bigfloat& ly, bigfloat& lz, bigfloat& d) const
{
	bigfloat p[4], q[4];
	ip.getBigfloatLambda(p[0], p[1], p[2], p[3]);
	iq.getBigfloatLambda(q[0], q[1], q[2], q[3]);
	lambda3d_LNC<bigfloat, bigfloat, bigfloat>(p, q, T(), lx, ly, lz, d);
}

inline void implicitPoint3D_BPT::getBigfloatLambda(bigfloat& lx, bigfloat& ly, bigfloat& lz, bigfloat& d) const
{
	bigfloat p[4], q[4], r[4];
	ip.getBigfloatLambda(p[0], p[1], p[2], p[3]);
	iq.getBigfloatLambda(q[0], q[1], q[2], q[3]);
	ir.getBigfloatLambda(r[0], r[1], r[2], r[3]);
	lambda3d_BPT<bigfloat, bigfloat, bigfloat>(p, q, r, U(), V(), lx, ly, lz, d);
}

inline void implicitPoint3D_TBC::getBigfloatLambda(bigfloat& lx, bigfloat& ly, bigfloat& lz, bigfloat& d) const
{
	bigfloat p[4], q[4], r[4], s[4];
	ip.getBigfloatLambda(p[0], p[1], p[2], p[3]);
	iq.getBigfloatLambda(q[0], q[1], q[2], q[3]);
	ir.getBigfloatLambda(r[0], r[1], r[2], r[3]);
	is.getBigfloatLambda(s[0], s[1], s[2], s[3]);
	lambda3d_TBC<bigfloat, bigfloat>(p, q, r, s, lx, ly, lz, d);
}


inline bool genericPoint::apapExplicit(explicitPoint2D& e) const
{
	if (isExplicit2D()) e = toExplicit2D();
	else {
		bigfloat l1x, l1y, d1;
		getBigfloatLambda(l1x, l1y, d1);
		const double lambda_x = l1x.get_d();
		const double lambda_y = l1y.get_d();
		const double lambda_d = d1.get_d();
		if (lambda_d == 0) return false;
		e = explicitPoint2D(lambda_x / lambda_d, lambda_y / lambda_d);
	}
	return true;
}

inline bool genericPoint::approxExplicit(explicitPoint2D& e) const
{
	if (isExplicit2D()) e = toExplicit2D();
	else {
		double lambda_x, lambda_y, lambda_d;
		interval_number ilx, ily, id;
		if (!getIntervalLambda(ilx, ily, id)) return apapExplicit(e);
		else
		{
			lambda_x = ilx.sup() + ilx.inf();
			lambda_y = ily.sup() + ily.inf();
			lambda_d = id.sup() + id.inf();
		}
		e = explicitPoint2D(lambda_x / lambda_d, lambda_y / lambda_d);
	}
	return true;
}

inline bool genericPoint::apapExplicit(explicitPoint3D& e) const
{
	if (isExplicit3D()) e = toExplicit3D();
	else {
		bigfloat l1z, l1x, l1y, d1;
		getBigfloatLambda(l1x, l1y, l1z, d1);
		const double lambda_x = l1x.get_d();
		const double lambda_y = l1y.get_d();
		const double lambda_z = l1z.get_d();
		const double lambda_d = d1.get_d();
		if (lambda_d == 0) return false;
		e = explicitPoint3D(lambda_x / lambda_d, lambda_y / lambda_d, lambda_z / lambda_d);
	}
	return true;
}


inline bool genericPoint::approxExplicit(explicitPoint3D& e) const
{
	if (isExplicit3D()) e = toExplicit3D();
	else {
		double lambda_x, lambda_y, lambda_z, lambda_d;
		interval_number ilx, ily, ilz, id;
		if (!getIntervalLambda(ilx, ily, ilz, id)) return apapExplicit(e);
		else
		{
			lambda_x = ilx.sup() + ilx.inf();
			lambda_y = ily.sup() + ily.inf();
			lambda_z = ilz.sup() + ilz.inf();
			lambda_d = id.sup() + id.inf();
		}
		e = explicitPoint3D(lambda_x / lambda_d, lambda_y / lambda_d, lambda_z / lambda_d);
	}
	return true;
}


inline bool genericPoint::getApproxXYCoordinates(double& x, double& y, bool apap) const
{
	if (is2D())
	{
		explicitPoint2D op;
		if (apap && !apapExplicit(op)) return false;
		if (!apap && !approxExplicit(op)) return false;
		x = op.X(); y = op.Y();
		return true;
	}
	if (is3D())
	{
		explicitPoint3D op;
		if (apap && !apapExplicit(op)) return false;
		if (!apap && !approxExplicit(op)) return false;
		x = op.X(); y = op.Y();
		return true;
	}
	ip_error("genericPoint::getApproxXYCoordinates - should not happen\n");
	return false;
}

inline bool genericPoint::getApproxXYZCoordinates(double& x, double& y, double& z, bool apap) const
{
	if (is3D())
	{
		explicitPoint3D op;
		if (apap && !apapExplicit(op)) return false;
		if (!apap && !approxExplicit(op)) return false;
		x = op.X(); y = op.Y(); z = op.Z();
		return true;
	}
	ip_error("genericPoint::getApproxXYZCoordinates - should not happen\n");
	return false;
}

inline bool genericPoint::getExactXYCoordinates(bigrational& x, bigrational& y) const
{
	if (isExplicit2D()) return toExplicit2D().getExactXYCoordinates(x, y);
	else if (isSSI()) return toSSI().getExactXYCoordinates(x, y);
	else ip_error("genericPoint::getExactXYCoordinates - should not happen\n");
	return false;
}

inline bool genericPoint::getExactXYZCoordinates(bigrational& x, bigrational& y, bigrational& z) const
{
	if (isExplicit3D()) return toExplicit3D().getExactXYZCoordinates(x, y, z);
	else if (isLPI()) return toLPI().getExactXYZCoordinates(x, y, z);
	else if (isTPI()) return toTPI().getExactXYZCoordinates(x, y, z);
	else if (isLNC()) return toLNC().getExactXYZCoordinates(x, y, z);
	else if (isBPT()) return toBPT().getExactXYZCoordinates(x, y, z);
	else if (isTBC()) return toTBC().getExactXYZCoordinates(x, y, z);
	else if (isExplicit2D()) { z = bigfloat(0); return toExplicit2D().getExactXYCoordinates(x, y); }
	else if (isSSI()) { z = bigfloat(0); return toSSI().getExactXYCoordinates(x, y); }
	else ip_error("genericPoint::getExactXYZCoordinates - should not happen\n");
	return false;
}

inline bool implicitPoint2D_SSI::getExactXYCoordinates(bigrational& x, bigrational& y) const
{
	bigfloat lx, ly, d;
	getBigfloatLambda(lx, ly, d);
	if (sgn(d) == 0) return false;
	const bigrational rd(d);
	x = bigrational(lx) / rd;
	y = bigrational(ly) / rd;
	return true;
}

inline bool implicitPoint3D_LPI::getExactXYZCoordinates(bigrational& x, bigrational& y, bigrational& z) const
{
	bigfloat lx, ly, lz, d;
	getBigfloatLambda(lx, ly, lz, d);
	if (sgn(d) == 0) return false;
	const bigrational rd(d);
	x = bigrational(lx) / rd;
	y = bigrational(ly) / rd;
	z = bigrational(lz) / rd;
	return true;
}

inline bool implicitPoint3D_TPI::getExactXYZCoordinates(bigrational& x, bigrational& y, bigrational& z) const
{
	bigfloat lx, ly, lz, d;
	getBigfloatLambda(lx, ly, lz, d);
	if (sgn(d) == 0) return false;
	const bigrational rd(d);
	x = bigrational(lx) / rd;
	y = bigrational(ly) / rd;
	z = bigrational(lz) / rd;
	return true;
}

inline bool implicitPoint3D_LNC::getExactXYZCoordinates(bigrational& x, bigrational& y, bigrational& z) const
{
	bigfloat lx, ly, lz, d;
	getBigfloatLambda(lx, ly, lz, d);
	x = bigrational(lx);
	y = bigrational(ly);
	z = bigrational(lz);
	return true;
}

inline bool implicitPoint3D_BPT::getExactXYZCoordinates(bigrational& x, bigrational& y, bigrational& z) const
{
	bigfloat lx, ly, lz, d;
	getBigfloatLambda(lx, ly, lz, d);
	x = bigrational(lx);
	y = bigrational(ly);
	z = bigrational(lz);
	return true;
}

inline bool implicitPoint3D_TBC::getExactXYZCoordinates(bigrational& x, bigrational& y, bigrational& z) const
{
	bigfloat lx, ly, lz, d;
	getBigfloatLambda(lx, ly, lz, d);
	x = bigrational(lx);
	y = bigrational(ly);
	z = bigrational(lz);
	return true;
}

inline ostream& operator<<(ostream& os, const genericPoint& p)
{
	if (p.isExplicit2D()) return os << p.toExplicit2D();
	else if (p.isExplicit3D()) return os << p.toExplicit3D();
	else if (p.isSSI()) return os << p.toSSI();
	else if (p.isLPI()) return os << p.toLPI();
	else if (p.isTPI()) return os << p.toTPI();
	else if (p.isLNC()) return os << p.toLNC();
	else if (p.isBPT()) return os << p.toBPT();
	else if (p.isTBC()) return os << p.toTBC();
	else ip_error("genericPoint::operator<< - should not happen\n");
	return os;
}

inline int maxComponentInTriangleNormal_filtered(double ov1x, double ov1y, double ov1z, double ov2x, double ov2y, double ov2z, double ov3x, double ov3y, double ov3z)
{
	double v3x = ov3x - ov2x;
	double v3y = ov3y - ov2y;
	double v3z = ov3z - ov2z;
	double v2x = ov2x - ov1x;
	double v2y = ov2y - ov1y;
	double v2z = ov2z - ov1z;
	double nvx1 = v2y * v3z;
	double nvx2 = v2z * v3y;
	double nvx = nvx1 - nvx2;
	double nvy1 = v3x * v2z;
	double nvy2 = v3z * v2x;
	double nvy = nvy1 - nvy2;
	double nvz1 = v2x * v3y;
	double nvz2 = v2y * v3x;
	double nvz = nvz1 - nvz2;

	double _tmp_fabs, max_var = 0;
	if ((_tmp_fabs = fabs(v3x)) > max_var) max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(v3y)) > max_var) max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(v3z)) > max_var) max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(v2x)) > max_var) max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(v2y)) > max_var) max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(v2z)) > max_var) max_var = _tmp_fabs;
	double epsilon = 8.88395e-016 * max_var * max_var;

	double nvxc = fabs(nvx);
	double nvyc = fabs(nvy);
	double nvzc = fabs(nvz);
	double nv = nvxc;
	if (nvyc > nv) nv = nvyc;
	if (nvzc > nv) nv = nvzc;

	if (nv > epsilon)
	{
		if (nv == nvxc) return 0;
		if (nv == nvyc) return 1;
		if (nv == nvzc) return 2;
	}
	return -1;
}

inline int maxComponentInTriangleNormal_exact(double ov1x, double ov1y, double ov1z, double ov2x, double ov2y, double ov2z, double ov3x, double ov3y, double ov3z)
{
	expansionObject o;
	double v3x[2];
	o.two_Diff(ov3x, ov2x, v3x);
	double v3y[2];
	o.two_Diff(ov3y, ov2y, v3y);
	double v3z[2];
	o.two_Diff(ov3z, ov2z, v3z);
	double v2x[2];
	o.two_Diff(ov2x, ov1x, v2x);
	double v2y[2];
	o.two_Diff(ov2y, ov1y, v2y);
	double v2z[2];
	o.two_Diff(ov2z, ov1z, v2z);
	double nvx1[8];
	o.Two_Two_Prod(v2y, v3z, nvx1);
	double nvx2[8];
	o.Two_Two_Prod(v2z, v3y, nvx2);
	double nvx[16];
	int nvx_len = o.Gen_Diff(8, nvx1, 8, nvx2, nvx);
	double nvy1[8];
	o.Two_Two_Prod(v3x, v2z, nvy1);
	double nvy2[8];
	o.Two_Two_Prod(v3z, v2x, nvy2);
	double nvy[16];
	int nvy_len = o.Gen_Diff(8, nvy1, 8, nvy2, nvy);
	double nvz1[8];
	o.Two_Two_Prod(v2x, v3y, nvz1);
	double nvz2[8];
	o.Two_Two_Prod(v2y, v3x, nvz2);
	double nvz[16];
	int nvz_len = o.Gen_Diff(8, nvz1, 8, nvz2, nvz);

	double nvxc = fabs(nvx[nvx_len - 1]);
	double nvyc = fabs(nvy[nvy_len - 1]);
	double nvzc = fabs(nvz[nvz_len - 1]);
	double nv = nvxc;
	if (nvyc > nv) nv = nvyc;
	if (nvzc > nv) return 2;
	if (nv == nvxc) return 0;
	return 1;
}

inline int genericPoint::maxComponentInTriangleNormal(double ov1x, double ov1y, double ov1z, double ov2x, double ov2y, double ov2z, double ov3x, double ov3y, double ov3z)
{
	int ret;
	if ((ret = maxComponentInTriangleNormal_filtered(ov1x, ov1y, ov1z, ov2x, ov2y, ov2z, ov3x, ov3y, ov3z)) >= 0) return ret;
	return maxComponentInTriangleNormal_exact(ov1x, ov1y, ov1z, ov2x, ov2y, ov2z, ov3x, ov3y, ov3z);
}


/////////////////////////////////////////////
//
// Derived predicates
//
/////////////////////////////////////////////

inline bool genericPoint::innerSegmentsCross(const genericPoint& A, const genericPoint& B, const genericPoint& P, const genericPoint& Q)
{
	int o11, o12, o21, o22;

	o11 = orient2Dxy(P, A, B);
	o12 = orient2Dxy(Q, B, A);
	o21 = orient2Dxy(A, P, Q);
	o22 = orient2Dxy(B, Q, P);
	if (o11 || o21 || o12 || o22) return (o11 == o12 && o21 == o22);

	o11 = orient2Dyz(P, A, B);
	o12 = orient2Dyz(Q, B, A);
	o21 = orient2Dyz(A, P, Q);
	o22 = orient2Dyz(B, Q, P);
	if (o11 || o21 || o12 || o22) return (o11 == o12 && o21 == o22);

	o11 = orient2Dzx(P, A, B);
	o12 = orient2Dzx(Q, B, A);
	o21 = orient2Dzx(A, P, Q);
	o22 = orient2Dzx(B, Q, P);
	if (o11 || o21 || o12 || o22) return (o11 == o12 && o21 == o22);

	return false;
}

inline bool genericPoint::segmentsCross(const genericPoint& A, const genericPoint& B, const genericPoint& P, const genericPoint& Q)
{
	int o11, o12, o21, o22;

	o11 = orient2Dxy(P, A, B);
	o12 = orient2Dxy(Q, B, A);
	o21 = orient2Dxy(A, P, Q);
	o22 = orient2Dxy(B, Q, P);
	if ((o11 || o12) && (o11 * o12 >= 0) && (o21 || o22) && (o21 * o22 >= 0)) return true;

	o11 = orient2Dyz(P, A, B);
	o12 = orient2Dyz(Q, B, A);
	o21 = orient2Dyz(A, P, Q);
	o22 = orient2Dyz(B, Q, P);
	if ((o11 || o12) && (o11 * o12 >= 0) && (o21 || o22) && (o21 * o22 >= 0)) return true;

	o11 = orient2Dzx(P, A, B);
	o12 = orient2Dzx(Q, B, A);
	o21 = orient2Dzx(A, P, Q);
	o22 = orient2Dzx(B, Q, P);
	if ((o11 || o12) && (o11 * o12 >= 0) && (o21 || o22) && (o21 * o22 >= 0)) return true;

	return false;
}

inline bool genericPoint::innerSegmentCrossesInnerTriangle(const genericPoint& s1, const genericPoint& s2, const genericPoint& v1, const genericPoint& v2, const genericPoint& v3)
{
	int o1 = orient3D(s1, v1, v2, v3); if (o1 == 0) return false;
	int o2 = orient3D(s2, v1, v2, v3); if (o2 == 0) return false;

	if ((o1 > 0 && o2 > 0) || (o1 < 0 && o2 < 0)) return false;
	o1 = orient3D(s1, s2, v1, v2);
	o2 = orient3D(s1, s2, v2, v3);
	if ((o1 >= 0 && o2 <= 0) || (o1 <= 0 && o2 >= 0)) return false;
	int o3 = orient3D(s1, s2, v3, v1);
	if ((o1 >= 0 && o3 <= 0) || (o1 <= 0 && o3 >= 0)) return false;
	if ((o2 >= 0 && o3 <= 0) || (o2 <= 0 && o3 >= 0)) return false;
	return true;
}

inline bool genericPoint::pointInInnerSegment(const genericPoint& p, const genericPoint& v1, const genericPoint& v2)
{
	if (misaligned(p, v1, v2)) return false;

	int lt2, lt3;
	lt2 = lessThanOnX(v1, p);
	lt3 = lessThanOnX(p, v2);
	if (lt2) return (lt2 == lt3);
	lt2 = lessThanOnY(v1, p);
	lt3 = lessThanOnY(p, v2);
	if (lt2) return (lt2 == lt3);
	lt2 = lessThanOnZ(v1, p);
	lt3 = lessThanOnZ(p, v2);
	if (lt2) return (lt2 == lt3);
	return false;
}

inline bool genericPoint::pointInSegment(const genericPoint& p, const genericPoint& v1, const genericPoint& v2)
{
	if (misaligned(p, v1, v2)) return false;

	int lt2x = lessThanOnX(v1, p);
	int lt3x = lessThanOnX(p, v2);
	if (lt2x && lt3x) return (lt2x == lt3x);
	int lt2y = lessThanOnY(v1, p);
	int lt3y = lessThanOnY(p, v2);
	if (lt2y && lt3y) return (lt2y == lt3y);
	int lt2z = lessThanOnZ(v1, p);
	int lt3z = lessThanOnZ(p, v2);
	if (lt2z && lt3z) return (lt2z == lt3z);

	return ((lt2x == 0 && lt2y == 0 && lt2z == 0) || (lt3x == 0 && lt3y == 0 && lt3z == 0));
}

inline bool genericPoint::pointInTriangle(const genericPoint& P, const genericPoint& A, const genericPoint& B, const genericPoint& C)
{
	int o1, o2, o3;
	o1 = orient2Dxy(P, A, B);
	o2 = orient2Dxy(P, B, C);
	o3 = orient2Dxy(P, C, A);
	if (o1 || o2 || o3) return ((o1 >= 0 && o2 >= 0 && o3 >= 0) || (o1 <= 0 && o2 <= 0 && o3 <= 0));
	o1 = orient2Dyz(P, A, B);
	o2 = orient2Dyz(P, B, C);
	o3 = orient2Dyz(P, C, A);
	if (o1 || o2 || o3) return ((o1 >= 0 && o2 >= 0 && o3 >= 0) || (o1 <= 0 && o2 <= 0 && o3 <= 0));
	o1 = orient2Dzx(P, A, B);
	o2 = orient2Dzx(P, B, C);
	o3 = orient2Dzx(P, C, A);
	return ((o1 >= 0 && o2 >= 0 && o3 >= 0) || (o1 <= 0 && o2 <= 0 && o3 <= 0));
}


inline bool genericPoint::pointInTriangle(const genericPoint& P, const genericPoint& A, const genericPoint& B, const genericPoint& C, int& o1, int& o2, int& o3)
{
	o1 = orient2Dxy(P, A, B);
	o2 = orient2Dxy(P, B, C);
	o3 = orient2Dxy(P, C, A);
	if (o1 || o2 || o3) return ((o1 >= 0 && o2 >= 0 && o3 >= 0) || (o1 <= 0 && o2 <= 0 && o3 <= 0));
	o1 = orient2Dyz(P, A, B);
	o2 = orient2Dyz(P, B, C);
	o3 = orient2Dyz(P, C, A);
	if (o1 || o2 || o3) return ((o1 >= 0 && o2 >= 0 && o3 >= 0) || (o1 <= 0 && o2 <= 0 && o3 <= 0));
	o1 = orient2Dzx(P, A, B);
	o2 = orient2Dzx(P, B, C);
	o3 = orient2Dzx(P, C, A);
	return ((o1 >= 0 && o2 >= 0 && o3 >= 0) || (o1 <= 0 && o2 <= 0 && o3 <= 0));
}

inline bool genericPoint::pointInInnerTriangle(const genericPoint& P, const genericPoint& A, const genericPoint& B, const genericPoint& C)
{
	int o1, o2, o3;
	o1 = orient2Dxy(P, A, B);
	o2 = orient2Dxy(P, B, C);
	o3 = orient2Dxy(P, C, A);
	if (o1 || o2 || o3) return ((o1 > 0 && o2 > 0 && o3 > 0) || (o1 < 0 && o2 < 0 && o3 < 0));
	o1 = orient2Dyz(P, A, B);
	o2 = orient2Dyz(P, B, C);
	o3 = orient2Dyz(P, C, A);
	if (o1 || o2 || o3) return ((o1 > 0 && o2 > 0 && o3 > 0) || (o1 < 0 && o2 < 0 && o3 < 0));
	o1 = orient2Dzx(P, A, B);
	o2 = orient2Dzx(P, B, C);
	o3 = orient2Dzx(P, C, A);
	return ((o1 > 0 && o2 > 0 && o3 > 0) || (o1 < 0 && o2 < 0 && o3 < 0));
}

inline bool genericPoint::lineCrossesInnerTriangle(const genericPoint& s1, const genericPoint& s2, const genericPoint& v1, const genericPoint& v2, const genericPoint& v3)
{
	const int o1 = genericPoint::orient3D(s1, s2, v1, v2);
	const int o2 = genericPoint::orient3D(s1, s2, v2, v3);
	if ((o1 >= 0 && o2 <= 0) || (o1 <= 0 && o2 >= 0)) return false;
	const int o3 = genericPoint::orient3D(s1, s2, v3, v1);
	if ((o1 >= 0 && o3 <= 0) || (o1 <= 0 && o3 >= 0)) return false;
	if ((o2 >= 0 && o3 <= 0) || (o2 <= 0 && o3 >= 0)) return false;
	return true;
}

inline bool genericPoint::lineCrossesTriangle(const genericPoint& s1, const genericPoint& s2, const genericPoint& v1, const genericPoint& v2, const genericPoint& v3)
{
	const int o1 = genericPoint::orient3D(s1, s2, v1, v2);
	const int o2 = genericPoint::orient3D(s1, s2, v2, v3);
	if ((o1 > 0 && o2 < 0) || (o1 < 0 && o2 > 0)) return false;
	const int o3 = genericPoint::orient3D(s1, s2, v3, v1);
	if ((o1 > 0 && o3 < 0) || (o1 < 0 && o3 > 0)) return false;
	if ((o2 > 0 && o3 < 0) || (o2 < 0 && o3 > 0)) return false;
	return true;
}

inline bool genericPoint::innerSegmentCrossesTriangle(const genericPoint& s1, const genericPoint& s2, const genericPoint& v1, const genericPoint& v2, const genericPoint& v3)
{
	const int o1 = genericPoint::orient3D(s1, v1, v2, v3); if (o1 == 0) return false;
	const int o2 = genericPoint::orient3D(s2, v1, v2, v3); if (o2 == 0) return false;

	if ((o1 > 0 && o2 > 0) || (o1 < 0 && o2 < 0)) return false;
	return lineCrossesTriangle(s1, s2, v1, v2, v3);
}



inline bool genericPoint::pointInInnerSegment(const genericPoint& p, const genericPoint& v1, const genericPoint& v2, int xyz)
{
	int lt2, lt3;
	if (xyz == 0)
	{
		if (orient2Dyz(p, v1, v2)) return false;
		lt2 = lessThanOnY(v1, p);
		lt3 = lessThanOnY(p, v2);
		if (lt2) return (lt2 == lt3);
		lt2 = lessThanOnZ(v1, p);
		lt3 = lessThanOnZ(p, v2);
	}
	else if (xyz == 1)
	{
		if (orient2Dzx(p, v1, v2)) return false;
		lt2 = lessThanOnX(v1, p);
		lt3 = lessThanOnX(p, v2);
		if (lt2) return (lt2 == lt3);
		lt2 = lessThanOnZ(v1, p);
		lt3 = lessThanOnZ(p, v2);
	}
	else
	{
		if (orient2Dxy(p, v1, v2)) return false;
		lt2 = lessThanOnX(v1, p);
		lt3 = lessThanOnX(p, v2);
		if (lt2) return (lt2 == lt3);
		lt2 = lessThanOnY(v1, p);
		lt3 = lessThanOnY(p, v2);
	}
	return (lt2 && (lt2 == lt3));
}

inline bool genericPoint::pointInSegment(const genericPoint& p, const genericPoint& v1, const genericPoint& v2, int xyz)
{
	int lt2, lt3, lt4, lt5;
	if (xyz == 0)
	{
		if (orient2Dyz(p, v1, v2)) return false;
		lt2 = lessThanOnY(v1, p);
		lt3 = lessThanOnY(p, v2);
		if (lt2 && lt3) return (lt2 == lt3);
		lt4 = lessThanOnZ(v1, p);
		lt5 = lessThanOnZ(p, v2);
	}
	else if (xyz == 1)
	{
		if (orient2Dzx(p, v1, v2)) return false;
		lt2 = lessThanOnX(v1, p);
		lt3 = lessThanOnX(p, v2);
		if (lt2 && lt3) return (lt2 == lt3);
		lt4 = lessThanOnZ(v1, p);
		lt5 = lessThanOnZ(p, v2);
	}
	else
	{
		if (orient2Dxy(p, v1, v2)) return false;
		lt2 = lessThanOnX(v1, p);
		lt3 = lessThanOnX(p, v2);
		if (lt2 && lt3) return (lt2 == lt3);
		lt4 = lessThanOnY(v1, p);
		lt5 = lessThanOnY(p, v2);
	}
	return ((lt2 == 0 && lt4 == 0) || (lt3 == 0 && lt5 == 0));
}

inline bool genericPoint::pointInInnerTriangle(const genericPoint& P, const genericPoint& A, const genericPoint& B, const genericPoint& C, int xyz)
{
	int o1, o2, o3;
	if (xyz == 2)
	{
		o1 = genericPoint::orient2Dxy(P, A, B);
		o2 = genericPoint::orient2Dxy(P, B, C);
		o3 = genericPoint::orient2Dxy(P, C, A);
	}
	else if (xyz == 0)
	{
		o1 = genericPoint::orient2Dyz(P, A, B);
		o2 = genericPoint::orient2Dyz(P, B, C);
		o3 = genericPoint::orient2Dyz(P, C, A);
	}
	else
	{
		o1 = genericPoint::orient2Dzx(P, A, B);
		o2 = genericPoint::orient2Dzx(P, B, C);
		o3 = genericPoint::orient2Dzx(P, C, A);
	}
	return ((o1 > 0 && o2 > 0 && o3 > 0) || (o1 < 0 && o2 < 0 && o3 < 0));
}

inline bool genericPoint::pointInTriangle(const genericPoint& P, const genericPoint& A, const genericPoint& B, const genericPoint& C, int xyz)
{
	int o1, o2, o3;
	if (xyz == 2)
	{
		o1 = genericPoint::orient2Dxy(P, A, B);
		o2 = genericPoint::orient2Dxy(P, B, C);
		o3 = genericPoint::orient2Dxy(P, C, A);
	}
	else if (xyz == 0)
	{
		o1 = genericPoint::orient2Dyz(P, A, B);
		o2 = genericPoint::orient2Dyz(P, B, C);
		o3 = genericPoint::orient2Dyz(P, C, A);
	}
	else
	{
		o1 = genericPoint::orient2Dzx(P, A, B);
		o2 = genericPoint::orient2Dzx(P, B, C);
		o3 = genericPoint::orient2Dzx(P, C, A);
	}
	return ((o1 >= 0 && o2 >= 0 && o3 >= 0) || (o1 <= 0 && o2 <= 0 && o3 <= 0));
}

inline bool genericPoint::innerSegmentsCross(const genericPoint& A, const genericPoint& B, const genericPoint& P, const genericPoint& Q, int xyz)
{
	int o11, o12, o21, o22;

	if (xyz == 2)
	{
		o11 = orient2Dxy(P, A, B);
		o12 = orient2Dxy(Q, B, A);
		o21 = orient2Dxy(A, P, Q);
		o22 = orient2Dxy(B, Q, P);
	}
	else if (xyz == 0)
	{
		o11 = orient2Dyz(P, A, B);
		o12 = orient2Dyz(Q, B, A);
		o21 = orient2Dyz(A, P, Q);
		o22 = orient2Dyz(B, Q, P);
	}
	else
	{
		o11 = orient2Dzx(P, A, B);
		o12 = orient2Dzx(Q, B, A);
		o21 = orient2Dzx(A, P, Q);
		o22 = orient2Dzx(B, Q, P);
	}

	return (o11 && o21 && o11 == o12 && o21 == o22);
}

inline bool genericPoint::segmentsCross(const genericPoint& A, const genericPoint& B, const genericPoint& P, const genericPoint& Q, int xyz)
{
	int o11, o12, o21, o22;

	if (xyz == 2)
	{
		o11 = orient2Dxy(P, A, B);
		o12 = orient2Dxy(Q, B, A);
		o21 = orient2Dxy(A, P, Q);
		o22 = orient2Dxy(B, Q, P);
	}
	else if (xyz == 0)
	{
		o11 = orient2Dyz(P, A, B);
		o12 = orient2Dyz(Q, B, A);
		o21 = orient2Dyz(A, P, Q);
		o22 = orient2Dyz(B, Q, P);
	}
	else
	{
		o11 = orient2Dzx(P, A, B);
		o12 = orient2Dzx(Q, B, A);
		o21 = orient2Dzx(A, P, Q);
		o22 = orient2Dzx(B, Q, P);
	}

	return ((o11 || o12) && (o11 * o12 >= 0) && (o21 || o22) && (o21 * o22 >= 0));
}
