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

/* This code was generated automatically. Do not edit unless you exactly   */
/* know what you are doing!                                                */

#pragma once

#include "numerics.h"

#pragma intrinsic(fabs)

static inline double ipow2(const double& d) { return d * d; }
static inline double ipow3(const double& d) { return d * ipow2(d); }

static inline interval_number ipow2(const interval_number& d) { return d.pow2(); }
static inline interval_number ipow3(const interval_number& d) { return d.pow3(); }

static inline expansion ipow2(const expansion& d) { return d.sqr(); }
static inline expansion ipow3(const expansion& d) { return d.sqr()*d; }

static inline bigfloat ipow2(const bigfloat& d) { return d * d; }
static inline bigfloat ipow3(const bigfloat& d) { return d * ipow2(d); }

inline int sgn(const interval_number& p) { return (p.isPositive()) ? (1) : ((p.isNegative()) ? (-1) : (0)); }


template<class PT, class T> static inline int dotProductSign2D_t(const PT& px, const PT& py, const PT& rx, const PT& ry, const PT& qx, const PT& qy) {
	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundUP();

	std::conditional_t<(std::is_same<expansion, T>::value), expansionPool, char> pool;

	if constexpr (std::is_same<expansion, T>::value) {
		expansion::initPool(&pool);
		feclearexcept(FE_ALL_EXCEPT);
	}

	const T lx = (px-qx);
	const T ly = (py-qy);
	const T gx = (rx-qx);
	const T gy = (ry-qy);
	const T d = ((lx*gx)+(ly*gy));

	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();

	if constexpr (std::is_same<expansion, T>::value) {
		if (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW)) return INT_MAX;
	}

	if constexpr (std::is_same<double, T>::value) {
		double _tmp_fabs;
		double max_var = 0.0;

		_tmp_fabs = fabs(lx); max_var = (_tmp_fabs > max_var)?(_tmp_fabs):(max_var);
		_tmp_fabs = fabs(ly); max_var = (_tmp_fabs > max_var)?(_tmp_fabs):(max_var);
		_tmp_fabs = fabs(gx); max_var = (_tmp_fabs > max_var)?(_tmp_fabs):(max_var);
		_tmp_fabs = fabs(gy); max_var = (_tmp_fabs > max_var)?(_tmp_fabs):(max_var);
		double epsilon = max_var;

		epsilon *= epsilon;
		epsilon *= 8.881784197001252e-16;

		return (d > epsilon) - (d < -epsilon);
	}
	else {
		return sgn(d);
	}
}

inline int dotProductSign2D(const double& px, const double& py, const double& rx, const double& ry, const double& qx, const double& qy) {
	int ret;
	if ((ret = dotProductSign2D_t<double, double>(px, py, rx, ry, qx, qy)) != 0) return ret;
	if ((ret = dotProductSign2D_t<s_expansion, expansion>(px, py, rx, ry, qx, qy)) != INT_MAX) return ret;
	return dotProductSign2D_t<bigfloat, bigfloat>(px, py, rx, ry, qx, qy);
}


template<class PT, class T> static inline int dotProductSign3D_t(const PT& px, const PT& py, const PT& pz, const PT& rx, const PT& ry, const PT& rz, const PT& qx, const PT& qy, const PT& qz) {
	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundUP();

	std::conditional_t<(std::is_same<expansion, T>::value), expansionPool, char> pool;

	if constexpr (std::is_same<expansion, T>::value) {
		expansion::initPool(&pool);
		feclearexcept(FE_ALL_EXCEPT);
	}

	const T lx = (px-qx);
	const T ly = (py-qy);
	const T lz = (pz-qz);
	const T gx = (rx-qx);
	const T gy = (ry-qy);
	const T gz = (rz-qz);
	const T d = ((lx*gx)+((ly*gy)+(lz*gz)));

	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();

	if constexpr (std::is_same<expansion, T>::value) {
		if (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW)) return INT_MAX;
	}

	if constexpr (std::is_same<double, T>::value) {
		double _tmp_fabs;
		double max_var = 0.0;

		_tmp_fabs = fabs(lx); max_var = (_tmp_fabs > max_var)?(_tmp_fabs):(max_var);
		_tmp_fabs = fabs(ly); max_var = (_tmp_fabs > max_var)?(_tmp_fabs):(max_var);
		_tmp_fabs = fabs(lz); max_var = (_tmp_fabs > max_var)?(_tmp_fabs):(max_var);
		_tmp_fabs = fabs(gx); max_var = (_tmp_fabs > max_var)?(_tmp_fabs):(max_var);
		_tmp_fabs = fabs(gy); max_var = (_tmp_fabs > max_var)?(_tmp_fabs):(max_var);
		_tmp_fabs = fabs(gz); max_var = (_tmp_fabs > max_var)?(_tmp_fabs):(max_var);
		double epsilon = max_var;

		epsilon *= epsilon;
		epsilon *= 1.443289932012704e-15;

		return (d > epsilon) - (d < -epsilon);
	}
	else {
		return sgn(d);
	}
}

inline int dotProductSign3D(const double& px, const double& py, const double& pz, const double& rx, const double& ry, const double& rz, const double& qx, const double& qy, const double& qz) {
	int ret;
	if ((ret = dotProductSign3D_t<double, double>(px, py, pz, rx, ry, rz, qx, qy, qz)) != 0) return ret;
	if ((ret = dotProductSign3D_t<s_expansion, expansion>(px, py, pz, rx, ry, rz, qx, qy, qz)) != INT_MAX) return ret;
	return dotProductSign3D_t<bigfloat, bigfloat>(px, py, pz, rx, ry, rz, qx, qy, qz);
}


template<class PT, class T> static inline int incircle_t(const PT& pax, const PT& pay, const PT& pbx, const PT& pby, const PT& pcx, const PT& pcy, const PT& pdx, const PT& pdy) {
	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundUP();

	std::conditional_t<(std::is_same<expansion, T>::value), expansionPool, char> pool;

	if constexpr (std::is_same<expansion, T>::value) {
		expansion::initPool(&pool);
		feclearexcept(FE_ALL_EXCEPT);
	}

	const T adx = (pax-pdx);
	const T ady = (pay-pdy);
	const T bdx = (pbx-pdx);
	const T bdy = (pby-pdy);
	const T cdx = (pcx-pdx);
	const T cdy = (pcy-pdy);
	const T abdet = ((adx*bdy)-(bdx*ady));
	const T bcdet = ((bdx*cdy)-(cdx*bdy));
	const T cadet = ((cdx*ady)-(adx*cdy));
	const T alift = (ipow2(adx)+ipow2(ady));
	const T blift = (ipow2(bdx)+ipow2(bdy));
	const T clift = (ipow2(cdx)+ipow2(cdy));
	const T L = ((alift*bcdet)+((blift*cadet)+(clift*abdet)));

	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();

	if constexpr (std::is_same<expansion, T>::value) {
		if (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW)) return INT_MAX;
	}

	if constexpr (std::is_same<double, T>::value) {
		double _tmp_fabs;
		double max_var = 0.0;

		_tmp_fabs = fabs(adx); max_var = (_tmp_fabs > max_var)?(_tmp_fabs):(max_var);
		_tmp_fabs = fabs(ady); max_var = (_tmp_fabs > max_var)?(_tmp_fabs):(max_var);
		_tmp_fabs = fabs(bdx); max_var = (_tmp_fabs > max_var)?(_tmp_fabs):(max_var);
		_tmp_fabs = fabs(bdy); max_var = (_tmp_fabs > max_var)?(_tmp_fabs):(max_var);
		_tmp_fabs = fabs(cdx); max_var = (_tmp_fabs > max_var)?(_tmp_fabs):(max_var);
		_tmp_fabs = fabs(cdy); max_var = (_tmp_fabs > max_var)?(_tmp_fabs):(max_var);
		double epsilon = max_var;

		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= 1.376676550535194e-14;

		return (L > epsilon) - (L < -epsilon);
	}
	else {
		return sgn(L);
	}
}

inline int incircle(const double& pax, const double& pay, const double& pbx, const double& pby, const double& pcx, const double& pcy, const double& pdx, const double& pdy) {
	int ret;
	if ((ret = incircle_t<double, double>(pax, pay, pbx, pby, pcx, pcy, pdx, pdy)) != 0) return ret;
	if ((ret = incircle_t<s_expansion, expansion>(pax, pay, pbx, pby, pcx, pcy, pdx, pdy)) != INT_MAX) return ret;
	return incircle_t<bigfloat, bigfloat>(pax, pay, pbx, pby, pcx, pcy, pdx, pdy);
}


template<class PT, class T> static inline int inGabrielSphere_t(const PT& qx, const PT& qy, const PT& qz, const PT& ax, const PT& ay, const PT& az, const PT& bx, const PT& by, const PT& bz, const PT& cx, const PT& cy, const PT& cz) {
	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundUP();

	std::conditional_t<(std::is_same<expansion, T>::value), expansionPool, char> pool;

	if constexpr (std::is_same<expansion, T>::value) {
		expansion::initPool(&pool);
		feclearexcept(FE_ALL_EXCEPT);
	}

	const T bax = (bx-ax);
	const T bay = (by-ay);
	const T baz = (bz-az);
	const T cax = (cx-ax);
	const T cay = (cy-ay);
	const T caz = (cz-az);
	const T qax = (qx-ax);
	const T qay = (qy-ay);
	const T qaz = (qz-az);
	const T ba2 = (ipow2(bax)+(ipow2(bay)+ipow2(baz)));
	const T ca2 = (ipow2(cax)+(ipow2(cay)+ipow2(caz)));
	const T abcx = ((cax*ba2)-(bax*ca2));
	const T abcy = ((cay*ba2)-(bay*ca2));
	const T abcz = ((caz*ba2)-(baz*ca2));
	const T crossbcx = ((bay*caz)-(baz*cay));
	const T crossbcy = ((baz*cax)-(bax*caz));
	const T crossbcz = ((bax*cay)-(bay*cax));
	const T ccax = ((abcy*crossbcz)-(abcz*crossbcy));
	const T ccay = ((abcz*crossbcx)-(abcx*crossbcz));
	const T ccaz = ((abcx*crossbcy)-(abcy*crossbcx));
	const T c2 = (ipow2(crossbcx)+(ipow2(crossbcy)+ipow2(crossbcz)));
	const T ret = ((qax*((qax*c2)-ccax))+((qay*((qay*c2)-ccay))+(qaz*((qaz*c2)-ccaz))));

	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();

	if constexpr (std::is_same<expansion, T>::value) {
		if (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW)) return INT_MAX;
	}

	if constexpr (std::is_same<double, T>::value) {
		double _tmp_fabs;
		double max_var = 0.0;

		_tmp_fabs = fabs(bax); max_var = (_tmp_fabs > max_var)?(_tmp_fabs):(max_var);
		_tmp_fabs = fabs(bay); max_var = (_tmp_fabs > max_var)?(_tmp_fabs):(max_var);
		_tmp_fabs = fabs(baz); max_var = (_tmp_fabs > max_var)?(_tmp_fabs):(max_var);
		_tmp_fabs = fabs(cax); max_var = (_tmp_fabs > max_var)?(_tmp_fabs):(max_var);
		_tmp_fabs = fabs(cay); max_var = (_tmp_fabs > max_var)?(_tmp_fabs):(max_var);
		_tmp_fabs = fabs(caz); max_var = (_tmp_fabs > max_var)?(_tmp_fabs):(max_var);
		_tmp_fabs = fabs(qax); max_var = (_tmp_fabs > max_var)?(_tmp_fabs):(max_var);
		_tmp_fabs = fabs(qay); max_var = (_tmp_fabs > max_var)?(_tmp_fabs):(max_var);
		_tmp_fabs = fabs(qaz); max_var = (_tmp_fabs > max_var)?(_tmp_fabs):(max_var);
		double epsilon = max_var;

		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 1.914024494453771e-13;

		return (ret > epsilon) - (ret < -epsilon);
	}
	else {
		return sgn(ret);
	}
}

inline int inGabrielSphere(const double& qx, const double& qy, const double& qz, const double& ax, const double& ay, const double& az, const double& bx, const double& by, const double& bz, const double& cx, const double& cy, const double& cz) {
	int ret;
	if ((ret = inGabrielSphere_t<double, double>(qx, qy, qz, ax, ay, az, bx, by, bz, cx, cy, cz)) != 0) return ret;
	return inGabrielSphere_t<bigfloat, bigfloat>(qx, qy, qz, ax, ay, az, bx, by, bz, cx, cy, cz);
}


template<class PT, class T> static inline int inSphere_t(const PT& pax, const PT& pay, const PT& paz, const PT& pbx, const PT& pby, const PT& pbz, const PT& pcx, const PT& pcy, const PT& pcz, const PT& pdx, const PT& pdy, const PT& pdz, const PT& pex, const PT& pey, const PT& pez) {
	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundUP();

	std::conditional_t<(std::is_same<expansion, T>::value), expansionPool, char> pool;

	if constexpr (std::is_same<expansion, T>::value) {
		expansion::initPool(&pool);
		feclearexcept(FE_ALL_EXCEPT);
	}

	const T aex = (pax-pex);
	const T aey = (pay-pey);
	const T aez = (paz-pez);
	const T alift = (ipow2(aex)+(ipow2(aey)+ipow2(aez)));
	const T bex = (pbx-pex);
	const T bey = (pby-pey);
	const T bez = (pbz-pez);
	const T blift = (ipow2(bex)+(ipow2(bey)+ipow2(bez)));
	const T cex = (pcx-pex);
	const T cey = (pcy-pey);
	const T cez = (pcz-pez);
	const T clift = (ipow2(cex)+(ipow2(cey)+ipow2(cez)));
	const T dex = (pdx-pex);
	const T dey = (pdy-pey);
	const T dez = (pdz-pez);
	const T dlift = (ipow2(dex)+(ipow2(dey)+ipow2(dez)));
	const T ab = ((aex*bey)-(bex*aey));
	const T bc = ((bex*cey)-(cex*bey));
	const T cd = ((cex*dey)-(dex*cey));
	const T da = ((dex*aey)-(aex*dey));
	const T ac = ((aex*cey)-(cex*aey));
	const T bd = ((bex*dey)-(dex*bey));
	const T abc = (((aez*bc)-(bez*ac))+(cez*ab));
	const T bcd = (((bez*cd)-(cez*bd))+(dez*bc));
	const T cda = ((cez*da)+((dez*ac)+(aez*cd)));
	const T dab = ((dez*ab)+((aez*bd)+(bez*da)));
	const T d = (((clift*dab)-(dlift*abc))+((alift*bcd)-(blift*cda)));

	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();

	if constexpr (std::is_same<expansion, T>::value) {
		if (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW)) return INT_MAX;
	}

	if constexpr (std::is_same<double, T>::value) {
		double _tmp_fabs;
		double max_var = 0.0;

		_tmp_fabs = fabs(aex); max_var = (_tmp_fabs > max_var)?(_tmp_fabs):(max_var);
		_tmp_fabs = fabs(aey); max_var = (_tmp_fabs > max_var)?(_tmp_fabs):(max_var);
		_tmp_fabs = fabs(aez); max_var = (_tmp_fabs > max_var)?(_tmp_fabs):(max_var);
		_tmp_fabs = fabs(bex); max_var = (_tmp_fabs > max_var)?(_tmp_fabs):(max_var);
		_tmp_fabs = fabs(bey); max_var = (_tmp_fabs > max_var)?(_tmp_fabs):(max_var);
		_tmp_fabs = fabs(bez); max_var = (_tmp_fabs > max_var)?(_tmp_fabs):(max_var);
		_tmp_fabs = fabs(cex); max_var = (_tmp_fabs > max_var)?(_tmp_fabs):(max_var);
		_tmp_fabs = fabs(cey); max_var = (_tmp_fabs > max_var)?(_tmp_fabs):(max_var);
		_tmp_fabs = fabs(cez); max_var = (_tmp_fabs > max_var)?(_tmp_fabs):(max_var);
		_tmp_fabs = fabs(dex); max_var = (_tmp_fabs > max_var)?(_tmp_fabs):(max_var);
		_tmp_fabs = fabs(dey); max_var = (_tmp_fabs > max_var)?(_tmp_fabs):(max_var);
		_tmp_fabs = fabs(dez); max_var = (_tmp_fabs > max_var)?(_tmp_fabs):(max_var);
		double epsilon = max_var;

		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= 1.145750161413162e-13;

		return (d > epsilon) - (d < -epsilon);
	}
	else {
		return sgn(d);
	}
}

inline int inSphere(const double& pax, const double& pay, const double& paz, const double& pbx, const double& pby, const double& pbz, const double& pcx, const double& pcy, const double& pcz, const double& pdx, const double& pdy, const double& pdz, const double& pex, const double& pey, const double& pez) {
	int ret;
	if ((ret = inSphere_t<double, double>(pax, pay, paz, pbx, pby, pbz, pcx, pcy, pcz, pdx, pdy, pdz, pex, pey, pez)) != 0) return ret;
	return inSphere_t<bigfloat, bigfloat>(pax, pay, paz, pbx, pby, pbz, pcx, pcy, pcz, pdx, pdy, pdz, pex, pey, pez);
}


template<class PT, class T> static inline int dotProductSign2D_EEI_t(const genericPoint& q, const PT& px, const PT& py, const PT& rx, const PT& ry) {
	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundUP();

	std::conditional_t<(std::is_same<expansion, T>::value), expansionPool, char> pool;

	if constexpr (std::is_same<expansion, T>::value) {
		expansion::initPool(&pool);
		feclearexcept(FE_ALL_EXCEPT);
	}

	T lqx, lqy, dq;
	if (!q.getLambda2D(lqx, lqy, dq)) {
		if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();
		return 0;
	}

	const T pxq = (px*dq);
	const T pyq = (py*dq);
	const T rxq = (rx*dq);
	const T ryq = (ry*dq);
	const T lx = (pxq-lqx);
	const T ly = (pyq-lqy);
	const T gx = (rxq-lqx);
	const T gy = (ryq-lqy);
	const T dx = (lx*gx);
	const T dy = (ly*gy);
	const T d = (dx+dy);

	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();

	if constexpr (std::is_same<expansion, T>::value) {
		if (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW)) return INT_MAX;
	}
	return sgn(d);
}

inline int dotProductSign2D_EEI(const genericPoint& q, const double& px, const double& py, const double& rx, const double& ry) {
	int ret;
	if ((ret = dotProductSign2D_EEI_t<interval_number, interval_number>(q, px, py, rx, ry)) != 0) return ret;
	return dotProductSign2D_EEI_t<bigfloat, bigfloat>(q, px, py, rx, ry);
}


template<class PT, class T> static inline int dotProductSign2D_IEE_t(const genericPoint& p, const PT& rx, const PT& ry, const PT& qx, const PT& qy) {
	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundUP();

	std::conditional_t<(std::is_same<expansion, T>::value), expansionPool, char> pool;

	if constexpr (std::is_same<expansion, T>::value) {
		expansion::initPool(&pool);
		feclearexcept(FE_ALL_EXCEPT);
	}

	T lpx, lpy, dp;
	if (!p.getLambda2D(lpx, lpy, dp)) {
		if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();
		return 0;
	}

	const T qxd = (qx*dp);
	const T qyd = (qy*dp);
	const T lx = (lpx-qxd);
	const T ly = (lpy-qyd);
	const T gx = (rx-qx);
	const T gy = (ry-qy);
	const T dx = (lx*gx);
	const T dy = (ly*gy);
	const T d = (dx+dy);

	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();

	if constexpr (std::is_same<expansion, T>::value) {
		if (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW)) return INT_MAX;
	}
	return sgn(d);
}

inline int dotProductSign2D_IEE(const genericPoint& p, const double& rx, const double& ry, const double& qx, const double& qy) {
	int ret;
	if ((ret = dotProductSign2D_IEE_t<interval_number, interval_number>(p, rx, ry, qx, qy)) != 0) return ret;
	return dotProductSign2D_IEE_t<bigfloat, bigfloat>(p, rx, ry, qx, qy);
}


template<class PT, class T> static inline int dotProductSign2D_IEI_t(const genericPoint& p, const genericPoint& q, const PT& rx, const PT& ry) {
	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundUP();

	std::conditional_t<(std::is_same<expansion, T>::value), expansionPool, char> pool;

	if constexpr (std::is_same<expansion, T>::value) {
		expansion::initPool(&pool);
		feclearexcept(FE_ALL_EXCEPT);
	}

	T lpx, lpy, dp;
	if (!p.getLambda2D(lpx, lpy, dp)) {
		if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();
		return 0;
	}

	T lqx, lqy, dq;
	if (!q.getLambda2D(lqx, lqy, dq)) {
		if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();
		return 0;
	}

	const T dqp = (dq*dp);
	const T pxq = (lpx*dqp);
	const T pyq = (lpy*dqp);
	const T rxq = (rx*dq);
	const T ryq = (ry*dq);
	const T lqxd = (lqx*dp);
	const T lqyd = (lqy*dp);
	const T lx = (pxq-lqxd);
	const T ly = (pyq-lqyd);
	const T gx = (rxq-lqx);
	const T gy = (ryq-lqy);
	const T dx = (lx*gx);
	const T dy = (ly*gy);
	const T d = (dx+dy);

	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();

	if constexpr (std::is_same<expansion, T>::value) {
		if (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW)) return INT_MAX;
	}
	return sgn(d);
}

inline int dotProductSign2D_IEI(const genericPoint& p, const genericPoint& q, const double& rx, const double& ry) {
	int ret;
	if ((ret = dotProductSign2D_IEI_t<interval_number, interval_number>(p, q, rx, ry)) != 0) return ret;
	return dotProductSign2D_IEI_t<bigfloat, bigfloat>(p, q, rx, ry);
}


template<class PT, class T> static inline int dotProductSign2D_IIE_t(const genericPoint& p, const genericPoint& r, const PT& qx, const PT& qy) {
	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundUP();

	std::conditional_t<(std::is_same<expansion, T>::value), expansionPool, char> pool;

	if constexpr (std::is_same<expansion, T>::value) {
		expansion::initPool(&pool);
		feclearexcept(FE_ALL_EXCEPT);
	}

	T lpx, lpy, dp;
	if (!p.getLambda2D(lpx, lpy, dp)) {
		if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();
		return 0;
	}

	T lrx, lry, dr;
	if (!r.getLambda2D(lrx, lry, dr)) {
		if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();
		return 0;
	}

	const T qxd = (qx*dp);
	const T qyd = (qy*dp);
	const T lx = (lpx-qxd);
	const T ly = (lpy-qyd);
	const T qxr = (qx*dr);
	const T qyr = (qy*dr);
	const T gx = (lrx-qxr);
	const T gy = (lry-qyr);
	const T dx = (lx*gx);
	const T dy = (ly*gy);
	const T d = (dx+dy);

	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();

	if constexpr (std::is_same<expansion, T>::value) {
		if (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW)) return INT_MAX;
	}
	return sgn(d);
}

inline int dotProductSign2D_IIE(const genericPoint& p, const genericPoint& r, const double& qx, const double& qy) {
	int ret;
	if ((ret = dotProductSign2D_IIE_t<interval_number, interval_number>(p, r, qx, qy)) != 0) return ret;
	return dotProductSign2D_IIE_t<bigfloat, bigfloat>(p, r, qx, qy);
}


template<class PT, class T> static inline int dotProductSign2D_III_t(const genericPoint& p, const genericPoint& r, const genericPoint& q) {
	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundUP();

	std::conditional_t<(std::is_same<expansion, T>::value), expansionPool, char> pool;

	if constexpr (std::is_same<expansion, T>::value) {
		expansion::initPool(&pool);
		feclearexcept(FE_ALL_EXCEPT);
	}

	T lpx, lpy, dp;
	if (!p.getLambda2D(lpx, lpy, dp)) {
		if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();
		return 0;
	}

	T lrx, lry, dr;
	if (!r.getLambda2D(lrx, lry, dr)) {
		if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();
		return 0;
	}

	T lqx, lqy, dq;
	if (!q.getLambda2D(lqx, lqy, dq)) {
		if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();
		return 0;
	}

	const T qxd = (lqx*dp);
	const T qyd = (lqy*dp);
	const T lpxq = (lpx*dq);
	const T lpyq = (lpy*dq);
	const T lx = (lpxq-qxd);
	const T ly = (lpyq-qyd);
	const T qxr = (lqx*dr);
	const T qyr = (lqy*dr);
	const T lrxq = (lrx*dq);
	const T lryq = (lry*dq);
	const T gx = (lrxq-qxr);
	const T gy = (lryq-qyr);
	const T dx = (lx*gx);
	const T dy = (ly*gy);
	const T d = (dx+dy);

	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();

	if constexpr (std::is_same<expansion, T>::value) {
		if (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW)) return INT_MAX;
	}
	return sgn(d);
}

inline int dotProductSign2D_III(const genericPoint& p, const genericPoint& r, const genericPoint& q) {
	int ret;
	if ((ret = dotProductSign2D_III_t<interval_number, interval_number>(p, r, q)) != 0) return ret;
	return dotProductSign2D_III_t<bigfloat, bigfloat>(p, r, q);
}


template<class PT, class T> static inline int dotProductSign3D_EEI_t(const genericPoint& q, const PT& px, const PT& py, const PT& pz, const PT& rx, const PT& ry, const PT& rz) {
	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundUP();

	std::conditional_t<(std::is_same<expansion, T>::value), expansionPool, char> pool;

	if constexpr (std::is_same<expansion, T>::value) {
		expansion::initPool(&pool);
		feclearexcept(FE_ALL_EXCEPT);
	}

	T lqx, lqy, lqz, dq;
	if (!q.getLambda3D(lqx, lqy, lqz, dq)) {
		if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();
		return 0;
	}

	const T pxq = (px*dq);
	const T pyq = (py*dq);
	const T pzq = (pz*dq);
	const T rxq = (rx*dq);
	const T ryq = (ry*dq);
	const T rzq = (rz*dq);
	const T lx = (pxq-lqx);
	const T ly = (pyq-lqy);
	const T lz = (pzq-lqz);
	const T gx = (rxq-lqx);
	const T gy = (ryq-lqy);
	const T gz = (rzq-lqz);
	const T dx = (lx*gx);
	const T dy = (ly*gy);
	const T dz = (lz*gz);
	const T d1 = (dx+dy);
	const T d = (d1+dz);

	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();

	if constexpr (std::is_same<expansion, T>::value) {
		if (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW)) return INT_MAX;
	}
	return sgn(d);
}

inline int dotProductSign3D_EEI(const genericPoint& q, const double& px, const double& py, const double& pz, const double& rx, const double& ry, const double& rz) {
	int ret;
	if ((ret = dotProductSign3D_EEI_t<interval_number, interval_number>(q, px, py, pz, rx, ry, rz)) != 0) return ret;
	return dotProductSign3D_EEI_t<bigfloat, bigfloat>(q, px, py, pz, rx, ry, rz);
}


template<class PT, class T> static inline int dotProductSign3D_IEE_t(const genericPoint& p, const PT& rx, const PT& ry, const PT& rz, const PT& qx, const PT& qy, const PT& qz) {
	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundUP();

	std::conditional_t<(std::is_same<expansion, T>::value), expansionPool, char> pool;

	if constexpr (std::is_same<expansion, T>::value) {
		expansion::initPool(&pool);
		feclearexcept(FE_ALL_EXCEPT);
	}

	T lpx, lpy, lpz, dp;
	if (!p.getLambda3D(lpx, lpy, lpz, dp)) {
		if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();
		return 0;
	}

	const T qxd = (qx*dp);
	const T qyd = (qy*dp);
	const T qzd = (qz*dp);
	const T lx = (lpx-qxd);
	const T ly = (lpy-qyd);
	const T lz = (lpz-qzd);
	const T gx = (rx-qx);
	const T gy = (ry-qy);
	const T gz = (rz-qz);
	const T dx = (lx*gx);
	const T dy = (ly*gy);
	const T dz = (lz*gz);
	const T d1 = (dx+dy);
	const T d = (d1+dz);

	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();

	if constexpr (std::is_same<expansion, T>::value) {
		if (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW)) return INT_MAX;
	}
	return sgn(d);
}

inline int dotProductSign3D_IEE(const genericPoint& p, const double& rx, const double& ry, const double& rz, const double& qx, const double& qy, const double& qz) {
	int ret;
	if ((ret = dotProductSign3D_IEE_t<interval_number, interval_number>(p, rx, ry, rz, qx, qy, qz)) != 0) return ret;
	return dotProductSign3D_IEE_t<bigfloat, bigfloat>(p, rx, ry, rz, qx, qy, qz);
}


template<class PT, class T> static inline int dotProductSign3D_IEI_t(const genericPoint& p, const genericPoint& q, const PT& rx, const PT& ry, const PT& rz) {
	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundUP();

	std::conditional_t<(std::is_same<expansion, T>::value), expansionPool, char> pool;

	if constexpr (std::is_same<expansion, T>::value) {
		expansion::initPool(&pool);
		feclearexcept(FE_ALL_EXCEPT);
	}

	T lpx, lpy, lpz, dp;
	if (!p.getLambda3D(lpx, lpy, lpz, dp)) {
		if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();
		return 0;
	}

	T lqx, lqy, lqz, dq;
	if (!q.getLambda3D(lqx, lqy, lqz, dq)) {
		if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();
		return 0;
	}

	const T dqp = (dq*dp);
	const T pxq = (lpx*dqp);
	const T pyq = (lpy*dqp);
	const T pzq = (lpz*dqp);
	const T rxq = (rx*dq);
	const T ryq = (ry*dq);
	const T rzq = (rz*dq);
	const T lqxd = (lqx*dp);
	const T lqyd = (lqy*dp);
	const T lqzd = (lqz*dp);
	const T lx = (pxq-lqxd);
	const T ly = (pyq-lqyd);
	const T lz = (pzq-lqzd);
	const T gx = (rxq-lqx);
	const T gy = (ryq-lqy);
	const T gz = (rzq-lqz);
	const T dx = (lx*gx);
	const T dy = (ly*gy);
	const T dz = (lz*gz);
	const T d1 = (dx+dy);
	const T d = (d1+dz);

	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();

	if constexpr (std::is_same<expansion, T>::value) {
		if (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW)) return INT_MAX;
	}
	return sgn(d);
}

inline int dotProductSign3D_IEI(const genericPoint& p, const genericPoint& q, const double& rx, const double& ry, const double& rz) {
	int ret;
	if ((ret = dotProductSign3D_IEI_t<interval_number, interval_number>(p, q, rx, ry, rz)) != 0) return ret;
	return dotProductSign3D_IEI_t<bigfloat, bigfloat>(p, q, rx, ry, rz);
}


template<class PT, class T> static inline int dotProductSign3D_IIE_t(const genericPoint& p, const genericPoint& r, const PT& qx, const PT& qy, const PT& qz) {
	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundUP();

	std::conditional_t<(std::is_same<expansion, T>::value), expansionPool, char> pool;

	if constexpr (std::is_same<expansion, T>::value) {
		expansion::initPool(&pool);
		feclearexcept(FE_ALL_EXCEPT);
	}

	T lpx, lpy, lpz, dp;
	if (!p.getLambda3D(lpx, lpy, lpz, dp)) {
		if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();
		return 0;
	}

	T lrx, lry, lrz, dr;
	if (!r.getLambda3D(lrx, lry, lrz, dr)) {
		if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();
		return 0;
	}

	const T qxd = (qx*dp);
	const T qyd = (qy*dp);
	const T qzd = (qz*dp);
	const T lx = (lpx-qxd);
	const T ly = (lpy-qyd);
	const T lz = (lpz-qzd);
	const T qxr = (qx*dr);
	const T qyr = (qy*dr);
	const T qzr = (qz*dr);
	const T gx = (lrx-qxr);
	const T gy = (lry-qyr);
	const T gz = (lrz-qzr);
	const T dx = (lx*gx);
	const T dy = (ly*gy);
	const T dz = (lz*gz);
	const T d1 = (dx+dy);
	const T d = (d1+dz);

	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();

	if constexpr (std::is_same<expansion, T>::value) {
		if (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW)) return INT_MAX;
	}
	return sgn(d);
}

inline int dotProductSign3D_IIE(const genericPoint& p, const genericPoint& r, const double& qx, const double& qy, const double& qz) {
	int ret;
	if ((ret = dotProductSign3D_IIE_t<interval_number, interval_number>(p, r, qx, qy, qz)) != 0) return ret;
	return dotProductSign3D_IIE_t<bigfloat, bigfloat>(p, r, qx, qy, qz);
}


template<class PT, class T> static inline int dotProductSign3D_III_t(const genericPoint& p, const genericPoint& r, const genericPoint& q) {
	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundUP();

	std::conditional_t<(std::is_same<expansion, T>::value), expansionPool, char> pool;

	if constexpr (std::is_same<expansion, T>::value) {
		expansion::initPool(&pool);
		feclearexcept(FE_ALL_EXCEPT);
	}

	T lpx, lpy, lpz, dp;
	if (!p.getLambda3D(lpx, lpy, lpz, dp)) {
		if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();
		return 0;
	}

	T lrx, lry, lrz, dr;
	if (!r.getLambda3D(lrx, lry, lrz, dr)) {
		if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();
		return 0;
	}

	T lqx, lqy, lqz, dq;
	if (!q.getLambda3D(lqx, lqy, lqz, dq)) {
		if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();
		return 0;
	}

	const T qxd = (lqx*dp);
	const T qyd = (lqy*dp);
	const T qzd = (lqz*dp);
	const T lpxq = (lpx*dq);
	const T lpyq = (lpy*dq);
	const T lpzq = (lpz*dq);
	const T lx = (lpxq-qxd);
	const T ly = (lpyq-qyd);
	const T lz = (lpzq-qzd);
	const T qxr = (lqx*dr);
	const T qyr = (lqy*dr);
	const T qzr = (lqz*dr);
	const T lrxq = (lrx*dq);
	const T lryq = (lry*dq);
	const T lrzq = (lrz*dq);
	const T gx = (lrxq-qxr);
	const T gy = (lryq-qyr);
	const T gz = (lrzq-qzr);
	const T dx = (lx*gx);
	const T dy = (ly*gy);
	const T dz = (lz*gz);
	const T d1 = (dx+dy);
	const T d = (d1+dz);

	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();

	if constexpr (std::is_same<expansion, T>::value) {
		if (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW)) return INT_MAX;
	}
	return sgn(d);
}

inline int dotProductSign3D_III(const genericPoint& p, const genericPoint& r, const genericPoint& q) {
	int ret;
	if ((ret = dotProductSign3D_III_t<interval_number, interval_number>(p, r, q)) != 0) return ret;
	return dotProductSign3D_III_t<bigfloat, bigfloat>(p, r, q);
}


template<class PT, class T> static inline int incirclexy_indirect_IEEE_t(const genericPoint& p1, const PT& pbx, const PT& pby, const PT& pcx, const PT& pcy, const PT& pdx, const PT& pdy) {
	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundUP();

	std::conditional_t<(std::is_same<expansion, T>::value), expansionPool, char> pool;

	if constexpr (std::is_same<expansion, T>::value) {
		expansion::initPool(&pool);
		feclearexcept(FE_ALL_EXCEPT);
	}

	T l1x, l1y, l1z, d1;
	if (!p1.getLambda3D(l1x, l1y, l1z, d1)) {
		if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();
		return 0;
	}

	const T pdxt = (pdx*d1);
	const T pdyt = (pdy*d1);
	const T adx = (l1x-pdxt);
	const T ady = (l1y-pdyt);
	const T bdx = (pbx-pdx);
	const T bdy = (pby-pdy);
	const T cdx = (pcx-pdx);
	const T cdy = (pcy-pdy);
	const T abdeta = (adx*bdy);
	const T abdetb = (bdx*ady);
	const T abdet = (abdeta-abdetb);
	const T bcdeta = (bdx*cdy);
	const T bcdetb = (cdx*bdy);
	const T bcdet = (bcdeta-bcdetb);
	const T cadeta = (cdx*ady);
	const T cadetb = (adx*cdy);
	const T cadet = (cadeta-cadetb);
	const T alifta = (adx*adx);
	const T aliftb = (ady*ady);
	const T alift = (alifta+aliftb);
	const T blifta = (bdx*bdx);
	const T bliftb = (bdy*bdy);
	const T blift = (blifta+bliftb);
	const T clifta = (cdx*cdx);
	const T cliftb = (cdy*cdy);
	const T clift = (clifta+cliftb);
	const T la = (alift*bcdet);
	const T lbt = (blift*cadet);
	const T lb = (lbt*d1);
	const T lct = (clift*abdet);
	const T lc = (lct*d1);
	const T lab = (la+lb);
	const T L = (lab+lc);

	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();

	if constexpr (std::is_same<expansion, T>::value) {
		if (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW)) return INT_MAX;
	}
	return sgn(L);
}

inline int incirclexy_indirect_IEEE(const genericPoint& p1, const double& pbx, const double& pby, const double& pcx, const double& pcy, const double& pdx, const double& pdy) {
	int ret;
	if ((ret = incirclexy_indirect_IEEE_t<interval_number, interval_number>(p1, pbx, pby, pcx, pcy, pdx, pdy)) != 0) return ret;
	return incirclexy_indirect_IEEE_t<bigfloat, bigfloat>(p1, pbx, pby, pcx, pcy, pdx, pdy);
}


template<class PT, class T> static inline int incirclexy_indirect_IIEE_t(const genericPoint& p1, const genericPoint& p2, const PT& pcx, const PT& pcy, const PT& pdx, const PT& pdy) {
	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundUP();

	std::conditional_t<(std::is_same<expansion, T>::value), expansionPool, char> pool;

	if constexpr (std::is_same<expansion, T>::value) {
		expansion::initPool(&pool);
		feclearexcept(FE_ALL_EXCEPT);
	}

	T l1x, l1y, l1z, d1;
	if (!p1.getLambda3D(l1x, l1y, l1z, d1)) {
		if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();
		return 0;
	}

	T l2x, l2y, l2z, d2;
	if (!p2.getLambda3D(l2x, l2y, l2z, d2)) {
		if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();
		return 0;
	}

	const T pdx1 = (pdx*d1);
	const T pdy1 = (pdy*d1);
	const T adx = (l1x-pdx1);
	const T ady = (l1y-pdy1);
	const T pdx2 = (pdx*d2);
	const T pdy2 = (pdy*d2);
	const T bdx = (l2x-pdx2);
	const T bdy = (l2y-pdy2);
	const T cdx = (pcx-pdx);
	const T cdy = (pcy-pdy);
	const T abdeta = (adx*bdy);
	const T abdetb = (bdx*ady);
	const T abdet = (abdeta-abdetb);
	const T bcdeta = (bdx*cdy);
	const T bcdetb = (cdx*bdy);
	const T bcdet = (bcdeta-bcdetb);
	const T cadeta = (cdx*ady);
	const T cadetb = (adx*cdy);
	const T cadet = (cadeta-cadetb);
	const T alifta = (adx*adx);
	const T aliftb = (ady*ady);
	const T aliftt = (alifta+aliftb);
	const T alift = (aliftt*d2);
	const T blifta = (bdx*bdx);
	const T bliftb = (bdy*bdy);
	const T blift = (blifta+bliftb);
	const T clifta = (cdx*cdx);
	const T cliftb = (cdy*cdy);
	const T cliftt = (clifta+cliftb);
	const T clift = (cliftt*d2);
	const T la = (alift*bcdet);
	const T lb = (blift*cadet);
	const T lc = (clift*abdet);
	const T lab = (lc+lb);
	const T lab2 = (lab*d1);
	const T L = (lab2+la);

	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();

	if constexpr (std::is_same<expansion, T>::value) {
		if (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW)) return INT_MAX;
	}
	return sgn(L);
}

inline int incirclexy_indirect_IIEE(const genericPoint& p1, const genericPoint& p2, const double& pcx, const double& pcy, const double& pdx, const double& pdy) {
	int ret;
	if ((ret = incirclexy_indirect_IIEE_t<interval_number, interval_number>(p1, p2, pcx, pcy, pdx, pdy)) != 0) return ret;
	return incirclexy_indirect_IIEE_t<bigfloat, bigfloat>(p1, p2, pcx, pcy, pdx, pdy);
}


template<class PT, class T> static inline int incirclexy_indirect_IIIE_t(const genericPoint& p1, const genericPoint& p2, const genericPoint& p3, const PT& pdx, const PT& pdy) {
	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundUP();

	std::conditional_t<(std::is_same<expansion, T>::value), expansionPool, char> pool;

	if constexpr (std::is_same<expansion, T>::value) {
		expansion::initPool(&pool);
		feclearexcept(FE_ALL_EXCEPT);
	}

	T l1x, l1y, l1z, d1;
	if (!p1.getLambda3D(l1x, l1y, l1z, d1)) {
		if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();
		return 0;
	}

	T l2x, l2y, l2z, d2;
	if (!p2.getLambda3D(l2x, l2y, l2z, d2)) {
		if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();
		return 0;
	}

	T l3x, l3y, l3z, d3;
	if (!p3.getLambda3D(l3x, l3y, l3z, d3)) {
		if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();
		return 0;
	}

	const T pdx1 = (pdx*d1);
	const T pdy1 = (pdy*d1);
	const T adx = (l1x-pdx1);
	const T ady = (l1y-pdy1);
	const T pdx2 = (pdx*d2);
	const T pdy2 = (pdy*d2);
	const T bdx = (l2x-pdx2);
	const T bdy = (l2y-pdy2);
	const T pdx3 = (pdx*d3);
	const T pdy3 = (pdy*d3);
	const T cdx = (l3x-pdx3);
	const T cdy = (l3y-pdy3);
	const T abdeta = (adx*bdy);
	const T abdetb = (bdx*ady);
	const T abdet = (abdeta-abdetb);
	const T bcdeta = (bdx*cdy);
	const T bcdetb = (cdx*bdy);
	const T bcdet = (bcdeta-bcdetb);
	const T cadeta = (cdx*ady);
	const T cadetb = (adx*cdy);
	const T cadet = (cadeta-cadetb);
	const T alifta = (adx*adx);
	const T aliftb = (ady*ady);
	const T aliftt = (alifta+aliftb);
	const T alift2 = (aliftt*d2);
	const T alift = (alift2*d3);
	const T blifta = (bdx*bdx);
	const T bliftb = (bdy*bdy);
	const T bliftt = (blifta+bliftb);
	const T blift = (bliftt*d3);
	const T clifta = (cdx*cdx);
	const T cliftb = (cdy*cdy);
	const T cliftt = (clifta+cliftb);
	const T clift = (cliftt*d2);
	const T la = (alift*bcdet);
	const T lb = (blift*cadet);
	const T lc = (clift*abdet);
	const T lab2 = (lc+lb);
	const T lab = (lab2*d1);
	const T L = (lab+la);

	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();

	if constexpr (std::is_same<expansion, T>::value) {
		if (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW)) return INT_MAX;
	}
	return sgn(L);
}

inline int incirclexy_indirect_IIIE(const genericPoint& p1, const genericPoint& p2, const genericPoint& p3, const double& pdx, const double& pdy) {
	int ret;
	if ((ret = incirclexy_indirect_IIIE_t<interval_number, interval_number>(p1, p2, p3, pdx, pdy)) != 0) return ret;
	return incirclexy_indirect_IIIE_t<bigfloat, bigfloat>(p1, p2, p3, pdx, pdy);
}


template<class PT, class T> static inline int incirclexy_indirect_IIII_t(const genericPoint& p1, const genericPoint& p2, const genericPoint& p3, const genericPoint& p4) {
	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundUP();

	std::conditional_t<(std::is_same<expansion, T>::value), expansionPool, char> pool;

	if constexpr (std::is_same<expansion, T>::value) {
		expansion::initPool(&pool);
		feclearexcept(FE_ALL_EXCEPT);
	}

	T l1x, l1y, l1z, d1;
	if (!p1.getLambda3D(l1x, l1y, l1z, d1)) {
		if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();
		return 0;
	}

	T l2x, l2y, l2z, d2;
	if (!p2.getLambda3D(l2x, l2y, l2z, d2)) {
		if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();
		return 0;
	}

	T l3x, l3y, l3z, d3;
	if (!p3.getLambda3D(l3x, l3y, l3z, d3)) {
		if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();
		return 0;
	}

	T l4x, l4y, l4z, d4;
	if (!p4.getLambda3D(l4x, l4y, l4z, d4)) {
		if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();
		return 0;
	}

	const T l1xt = (l1x*d4);
	const T l1yt = (l1y*d4);
	const T l2xt = (l2x*d4);
	const T l2yt = (l2y*d4);
	const T l3xt = (l3x*d4);
	const T l3yt = (l3y*d4);
	const T l4x1 = (l4x*d1);
	const T l4y1 = (l4y*d1);
	const T adx = (l1xt-l4x1);
	const T ady = (l1yt-l4y1);
	const T l4x2 = (l4x*d2);
	const T l4y2 = (l4y*d2);
	const T bdx = (l2xt-l4x2);
	const T bdy = (l2yt-l4y2);
	const T l4x3 = (l4x*d3);
	const T l4y3 = (l4y*d3);
	const T cdx = (l3xt-l4x3);
	const T cdy = (l3yt-l4y3);
	const T abdeta = (adx*bdy);
	const T abdetb = (bdx*ady);
	const T abdet = (abdeta-abdetb);
	const T bcdeta = (bdx*cdy);
	const T bcdetb = (cdx*bdy);
	const T bcdet = (bcdeta-bcdetb);
	const T cadeta = (cdx*ady);
	const T cadetb = (adx*cdy);
	const T cadet = (cadeta-cadetb);
	const T alifta = (adx*adx);
	const T aliftb = (ady*ady);
	const T aliftt = (alifta+aliftb);
	const T alift2 = (aliftt*d2);
	const T alift = (alift2*d3);
	const T blifta = (bdx*bdx);
	const T bliftb = (bdy*bdy);
	const T bliftt = (blifta+bliftb);
	const T blift = (bliftt*d3);
	const T clifta = (cdx*cdx);
	const T cliftb = (cdy*cdy);
	const T cliftt = (clifta+cliftb);
	const T clift = (cliftt*d2);
	const T la = (alift*bcdet);
	const T lb = (blift*cadet);
	const T lc = (clift*abdet);
	const T lab2 = (lc+lb);
	const T lab = (lab2*d1);
	const T L = (lab+la);

	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();

	if constexpr (std::is_same<expansion, T>::value) {
		if (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW)) return INT_MAX;
	}
	return sgn(L);
}

inline int incirclexy_indirect_IIII(const genericPoint& p1, const genericPoint& p2, const genericPoint& p3, const genericPoint& p4) {
	int ret;
	if ((ret = incirclexy_indirect_IIII_t<interval_number, interval_number>(p1, p2, p3, p4)) != 0) return ret;
	return incirclexy_indirect_IIII_t<bigfloat, bigfloat>(p1, p2, p3, p4);
}


template<class PT, class T> static inline int incircle_indirect_IEEE_t(const genericPoint& p1, const PT& pbx, const PT& pby, const PT& pcx, const PT& pcy, const PT& pdx, const PT& pdy) {
	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundUP();

	std::conditional_t<(std::is_same<expansion, T>::value), expansionPool, char> pool;

	if constexpr (std::is_same<expansion, T>::value) {
		expansion::initPool(&pool);
		feclearexcept(FE_ALL_EXCEPT);
	}

	T l1x, l1y, d1;
	if (!p1.getLambda2D(l1x, l1y, d1)) {
		if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();
		return 0;
	}

	const T pdxt = (pdx*d1);
	const T pdyt = (pdy*d1);
	const T adx = (l1x-pdxt);
	const T ady = (l1y-pdyt);
	const T bdx = (pbx-pdx);
	const T bdy = (pby-pdy);
	const T cdx = (pcx-pdx);
	const T cdy = (pcy-pdy);
	const T abdeta = (adx*bdy);
	const T abdetb = (bdx*ady);
	const T abdet = (abdeta-abdetb);
	const T bcdeta = (bdx*cdy);
	const T bcdetb = (cdx*bdy);
	const T bcdet = (bcdeta-bcdetb);
	const T cadeta = (cdx*ady);
	const T cadetb = (adx*cdy);
	const T cadet = (cadeta-cadetb);
	const T alifta = (adx*adx);
	const T aliftb = (ady*ady);
	const T alift = (alifta+aliftb);
	const T blifta = (bdx*bdx);
	const T bliftb = (bdy*bdy);
	const T blift = (blifta+bliftb);
	const T clifta = (cdx*cdx);
	const T cliftb = (cdy*cdy);
	const T clift = (clifta+cliftb);
	const T la = (alift*bcdet);
	const T lbt = (blift*cadet);
	const T lb = (lbt*d1);
	const T lct = (clift*abdet);
	const T lc = (lct*d1);
	const T lab = (la+lb);
	const T L = (lab+lc);

	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();

	if constexpr (std::is_same<expansion, T>::value) {
		if (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW)) return INT_MAX;
	}
	return sgn(L);
}

inline int incircle_indirect_IEEE(const genericPoint& p1, const double& pbx, const double& pby, const double& pcx, const double& pcy, const double& pdx, const double& pdy) {
	int ret;
	if ((ret = incircle_indirect_IEEE_t<interval_number, interval_number>(p1, pbx, pby, pcx, pcy, pdx, pdy)) != 0) return ret;
	return incircle_indirect_IEEE_t<bigfloat, bigfloat>(p1, pbx, pby, pcx, pcy, pdx, pdy);
}


template<class PT, class T> static inline int incircle_indirect_IIEE_t(const genericPoint& p1, const genericPoint& p2, const PT& pcx, const PT& pcy, const PT& pdx, const PT& pdy) {
	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundUP();

	std::conditional_t<(std::is_same<expansion, T>::value), expansionPool, char> pool;

	if constexpr (std::is_same<expansion, T>::value) {
		expansion::initPool(&pool);
		feclearexcept(FE_ALL_EXCEPT);
	}

	T l1x, l1y, d1;
	if (!p1.getLambda2D(l1x, l1y, d1)) {
		if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();
		return 0;
	}

	T l2x, l2y, d2;
	if (!p2.getLambda2D(l2x, l2y, d2)) {
		if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();
		return 0;
	}

	const T pdx1 = (pdx*d1);
	const T pdy1 = (pdy*d1);
	const T adx = (l1x-pdx1);
	const T ady = (l1y-pdy1);
	const T pdx2 = (pdx*d2);
	const T pdy2 = (pdy*d2);
	const T bdx = (l2x-pdx2);
	const T bdy = (l2y-pdy2);
	const T cdx = (pcx-pdx);
	const T cdy = (pcy-pdy);
	const T abdeta = (adx*bdy);
	const T abdetb = (bdx*ady);
	const T abdet = (abdeta-abdetb);
	const T bcdeta = (bdx*cdy);
	const T bcdetb = (cdx*bdy);
	const T bcdet = (bcdeta-bcdetb);
	const T cadeta = (cdx*ady);
	const T cadetb = (adx*cdy);
	const T cadet = (cadeta-cadetb);
	const T alifta = (adx*adx);
	const T aliftb = (ady*ady);
	const T aliftt = (alifta+aliftb);
	const T alift = (aliftt*d2);
	const T blifta = (bdx*bdx);
	const T bliftb = (bdy*bdy);
	const T blift = (blifta+bliftb);
	const T clifta = (cdx*cdx);
	const T cliftb = (cdy*cdy);
	const T cliftt = (clifta+cliftb);
	const T clift = (cliftt*d2);
	const T la = (alift*bcdet);
	const T lb = (blift*cadet);
	const T lc = (clift*abdet);
	const T lab = (lc+lb);
	const T lab2 = (lab*d1);
	const T L = (lab2+la);

	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();

	if constexpr (std::is_same<expansion, T>::value) {
		if (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW)) return INT_MAX;
	}
	return sgn(L);
}

inline int incircle_indirect_IIEE(const genericPoint& p1, const genericPoint& p2, const double& pcx, const double& pcy, const double& pdx, const double& pdy) {
	int ret;
	if ((ret = incircle_indirect_IIEE_t<interval_number, interval_number>(p1, p2, pcx, pcy, pdx, pdy)) != 0) return ret;
	return incircle_indirect_IIEE_t<bigfloat, bigfloat>(p1, p2, pcx, pcy, pdx, pdy);
}


template<class PT, class T> static inline int incircle_indirect_IIIE_t(const genericPoint& p1, const genericPoint& p2, const genericPoint& p3, const PT& pdx, const PT& pdy) {
	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundUP();

	std::conditional_t<(std::is_same<expansion, T>::value), expansionPool, char> pool;

	if constexpr (std::is_same<expansion, T>::value) {
		expansion::initPool(&pool);
		feclearexcept(FE_ALL_EXCEPT);
	}

	T l1x, l1y, d1;
	if (!p1.getLambda2D(l1x, l1y, d1)) {
		if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();
		return 0;
	}

	T l2x, l2y, d2;
	if (!p2.getLambda2D(l2x, l2y, d2)) {
		if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();
		return 0;
	}

	T l3x, l3y, d3;
	if (!p3.getLambda2D(l3x, l3y, d3)) {
		if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();
		return 0;
	}

	const T pdx1 = (pdx*d1);
	const T pdy1 = (pdy*d1);
	const T adx = (l1x-pdx1);
	const T ady = (l1y-pdy1);
	const T pdx2 = (pdx*d2);
	const T pdy2 = (pdy*d2);
	const T bdx = (l2x-pdx2);
	const T bdy = (l2y-pdy2);
	const T pdx3 = (pdx*d3);
	const T pdy3 = (pdy*d3);
	const T cdx = (l3x-pdx3);
	const T cdy = (l3y-pdy3);
	const T abdeta = (adx*bdy);
	const T abdetb = (bdx*ady);
	const T abdet = (abdeta-abdetb);
	const T bcdeta = (bdx*cdy);
	const T bcdetb = (cdx*bdy);
	const T bcdet = (bcdeta-bcdetb);
	const T cadeta = (cdx*ady);
	const T cadetb = (adx*cdy);
	const T cadet = (cadeta-cadetb);
	const T alifta = (adx*adx);
	const T aliftb = (ady*ady);
	const T aliftt = (alifta+aliftb);
	const T alift2 = (aliftt*d2);
	const T alift = (alift2*d3);
	const T blifta = (bdx*bdx);
	const T bliftb = (bdy*bdy);
	const T bliftt = (blifta+bliftb);
	const T blift = (bliftt*d3);
	const T clifta = (cdx*cdx);
	const T cliftb = (cdy*cdy);
	const T cliftt = (clifta+cliftb);
	const T clift = (cliftt*d2);
	const T la = (alift*bcdet);
	const T lb = (blift*cadet);
	const T lc = (clift*abdet);
	const T lab2 = (lc+lb);
	const T lab = (lab2*d1);
	const T L = (lab+la);

	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();

	if constexpr (std::is_same<expansion, T>::value) {
		if (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW)) return INT_MAX;
	}
	return sgn(L);
}

inline int incircle_indirect_IIIE(const genericPoint& p1, const genericPoint& p2, const genericPoint& p3, const double& pdx, const double& pdy) {
	int ret;
	if ((ret = incircle_indirect_IIIE_t<interval_number, interval_number>(p1, p2, p3, pdx, pdy)) != 0) return ret;
	return incircle_indirect_IIIE_t<bigfloat, bigfloat>(p1, p2, p3, pdx, pdy);
}


template<class PT, class T> static inline int incircle_indirect_IIII_t(const genericPoint& p1, const genericPoint& p2, const genericPoint& p3, const genericPoint& p4) {
	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundUP();

	std::conditional_t<(std::is_same<expansion, T>::value), expansionPool, char> pool;

	if constexpr (std::is_same<expansion, T>::value) {
		expansion::initPool(&pool);
		feclearexcept(FE_ALL_EXCEPT);
	}

	T l1x, l1y, d1;
	if (!p1.getLambda2D(l1x, l1y, d1)) {
		if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();
		return 0;
	}

	T l2x, l2y, d2;
	if (!p2.getLambda2D(l2x, l2y, d2)) {
		if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();
		return 0;
	}

	T l3x, l3y, d3;
	if (!p3.getLambda2D(l3x, l3y, d3)) {
		if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();
		return 0;
	}

	T l4x, l4y, d4;
	if (!p4.getLambda2D(l4x, l4y, d4)) {
		if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();
		return 0;
	}

	const T l1xt = (l1x*d4);
	const T l1yt = (l1y*d4);
	const T l2xt = (l2x*d4);
	const T l2yt = (l2y*d4);
	const T l3xt = (l3x*d4);
	const T l3yt = (l3y*d4);
	const T l4x1 = (l4x*d1);
	const T l4y1 = (l4y*d1);
	const T adx = (l1xt-l4x1);
	const T ady = (l1yt-l4y1);
	const T l4x2 = (l4x*d2);
	const T l4y2 = (l4y*d2);
	const T bdx = (l2xt-l4x2);
	const T bdy = (l2yt-l4y2);
	const T l4x3 = (l4x*d3);
	const T l4y3 = (l4y*d3);
	const T cdx = (l3xt-l4x3);
	const T cdy = (l3yt-l4y3);
	const T abdeta = (adx*bdy);
	const T abdetb = (bdx*ady);
	const T abdet = (abdeta-abdetb);
	const T bcdeta = (bdx*cdy);
	const T bcdetb = (cdx*bdy);
	const T bcdet = (bcdeta-bcdetb);
	const T cadeta = (cdx*ady);
	const T cadetb = (adx*cdy);
	const T cadet = (cadeta-cadetb);
	const T alifta = (adx*adx);
	const T aliftb = (ady*ady);
	const T aliftt = (alifta+aliftb);
	const T alift2 = (aliftt*d2);
	const T alift = (alift2*d3);
	const T blifta = (bdx*bdx);
	const T bliftb = (bdy*bdy);
	const T bliftt = (blifta+bliftb);
	const T blift = (bliftt*d3);
	const T clifta = (cdx*cdx);
	const T cliftb = (cdy*cdy);
	const T cliftt = (clifta+cliftb);
	const T clift = (cliftt*d2);
	const T la = (alift*bcdet);
	const T lb = (blift*cadet);
	const T lc = (clift*abdet);
	const T lab2 = (lc+lb);
	const T lab = (lab2*d1);
	const T L = (lab+la);

	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();

	if constexpr (std::is_same<expansion, T>::value) {
		if (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW)) return INT_MAX;
	}
	return sgn(L);
}

inline int incircle_indirect_IIII(const genericPoint& p1, const genericPoint& p2, const genericPoint& p3, const genericPoint& p4) {
	int ret;
	if ((ret = incircle_indirect_IIII_t<interval_number, interval_number>(p1, p2, p3, p4)) != 0) return ret;
	return incircle_indirect_IIII_t<bigfloat, bigfloat>(p1, p2, p3, p4);
}


template<class PT, class T> static inline int inGabrielSphere_EIEE_t(const genericPoint& a, const PT& qx, const PT& qy, const PT& qz, const PT& bx, const PT& by, const PT& bz, const PT& cx, const PT& cy, const PT& cz) {
	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundUP();

	std::conditional_t<(std::is_same<expansion, T>::value), expansionPool, char> pool;

	if constexpr (std::is_same<expansion, T>::value) {
		expansion::initPool(&pool);
		feclearexcept(FE_ALL_EXCEPT);
	}

	T l2x, l2y, l2z, d2;
	if (!a.getLambda3D(l2x, l2y, l2z, d2)) {
		if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();
		return 0;
	}

	const T bax = ((bx*d2)-l2x);
	const T bay = ((by*d2)-l2y);
	const T baz = ((bz*d2)-l2z);
	const T cax = ((cx*d2)-l2x);
	const T cay = ((cy*d2)-l2y);
	const T caz = ((cz*d2)-l2z);
	const T qax = ((qx*d2)-l2x);
	const T qay = ((qy*d2)-l2y);
	const T qaz = ((qz*d2)-l2z);
	const T ba2 = (ipow2(bax)+(ipow2(bay)+ipow2(baz)));
	const T ca2 = (ipow2(cax)+(ipow2(cay)+ipow2(caz)));
	const T abcx = ((cax*ba2)-(bax*ca2));
	const T abcy = ((cay*ba2)-(bay*ca2));
	const T abcz = ((caz*ba2)-(baz*ca2));
	const T crossbcx = ((bay*caz)-(baz*cay));
	const T crossbcy = ((baz*cax)-(bax*caz));
	const T crossbcz = ((bax*cay)-(bay*cax));
	const T ccax = ((abcy*crossbcz)-(abcz*crossbcy));
	const T ccay = ((abcz*crossbcx)-(abcx*crossbcz));
	const T ccaz = ((abcx*crossbcy)-(abcy*crossbcx));
	const T c2 = (ipow2(crossbcx)+(ipow2(crossbcy)+ipow2(crossbcz)));
	const T ret = ((qax*((qax*c2)-ccax))+((qay*((qay*c2)-ccay))+(qaz*((qaz*c2)-ccaz))));

	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();

	if constexpr (std::is_same<expansion, T>::value) {
		if (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW)) return INT_MAX;
	}
	return sgn(ret);
}

inline int inGabrielSphere_EIEE(const genericPoint& a, const double& qx, const double& qy, const double& qz, const double& bx, const double& by, const double& bz, const double& cx, const double& cy, const double& cz) {
	int ret;
	if ((ret = inGabrielSphere_EIEE_t<interval_number, interval_number>(a, qx, qy, qz, bx, by, bz, cx, cy, cz)) != 0) return ret;
	return inGabrielSphere_EIEE_t<bigfloat, bigfloat>(a, qx, qy, qz, bx, by, bz, cx, cy, cz);
}


template<class PT, class T> static inline int inGabrielSphere_EIIE_t(const genericPoint& a, const genericPoint& b, const PT& qx, const PT& qy, const PT& qz, const PT& cx, const PT& cy, const PT& cz) {
	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundUP();

	std::conditional_t<(std::is_same<expansion, T>::value), expansionPool, char> pool;

	if constexpr (std::is_same<expansion, T>::value) {
		expansion::initPool(&pool);
		feclearexcept(FE_ALL_EXCEPT);
	}

	T l2x, l2y, l2z, d2;
	if (!a.getLambda3D(l2x, l2y, l2z, d2)) {
		if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();
		return 0;
	}

	T l3x, l3y, l3z, d3;
	if (!b.getLambda3D(l3x, l3y, l3z, d3)) {
		if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();
		return 0;
	}

	const T bax = ((l3x*d2)-(l2x*d3));
	const T bay = ((l3y*d2)-(l2y*d3));
	const T baz = ((l3z*d2)-(l2z*d3));
	const T cax = ((cx*d2)-l2x);
	const T cay = ((cy*d2)-l2y);
	const T caz = ((cz*d2)-l2z);
	const T qax = ((qx*d2)-l2x);
	const T qay = ((qy*d2)-l2y);
	const T qaz = ((qz*d2)-l2z);
	const T ba2 = (ipow2(bax)+(ipow2(bay)+ipow2(baz)));
	const T ca2 = ((ipow2(cax)+(ipow2(cay)+ipow2(caz)))*d3);
	const T abcx = ((cax*ba2)-(bax*ca2));
	const T abcy = ((cay*ba2)-(bay*ca2));
	const T abcz = ((caz*ba2)-(baz*ca2));
	const T crossbcx = ((bay*caz)-(baz*cay));
	const T crossbcy = ((baz*cax)-(bax*caz));
	const T crossbcz = ((bax*cay)-(bay*cax));
	const T ccax = ((abcy*crossbcz)-(abcz*crossbcy));
	const T ccay = ((abcz*crossbcx)-(abcx*crossbcz));
	const T ccaz = ((abcx*crossbcy)-(abcy*crossbcx));
	const T c2 = ((ipow2(crossbcx)+(ipow2(crossbcy)+ipow2(crossbcz)))*d3);
	const T ret = ((qax*((qax*c2)-ccax))+((qay*((qay*c2)-ccay))+(qaz*((qaz*c2)-ccaz))));

	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();

	if constexpr (std::is_same<expansion, T>::value) {
		if (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW)) return INT_MAX;
	}
	return sgn(ret);
}

inline int inGabrielSphere_EIIE(const genericPoint& a, const genericPoint& b, const double& qx, const double& qy, const double& qz, const double& cx, const double& cy, const double& cz) {
	int ret;
	if ((ret = inGabrielSphere_EIIE_t<interval_number, interval_number>(a, b, qx, qy, qz, cx, cy, cz)) != 0) return ret;
	return inGabrielSphere_EIIE_t<bigfloat, bigfloat>(a, b, qx, qy, qz, cx, cy, cz);
}


template<class PT, class T> static inline int inGabrielSphere_EIII_t(const genericPoint& a, const genericPoint& b, const genericPoint& c, const PT& qx, const PT& qy, const PT& qz) {
	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundUP();

	std::conditional_t<(std::is_same<expansion, T>::value), expansionPool, char> pool;

	if constexpr (std::is_same<expansion, T>::value) {
		expansion::initPool(&pool);
		feclearexcept(FE_ALL_EXCEPT);
	}

	T l2x, l2y, l2z, d2;
	if (!a.getLambda3D(l2x, l2y, l2z, d2)) {
		if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();
		return 0;
	}

	T l3x, l3y, l3z, d3;
	if (!b.getLambda3D(l3x, l3y, l3z, d3)) {
		if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();
		return 0;
	}

	T l4x, l4y, l4z, d4;
	if (!c.getLambda3D(l4x, l4y, l4z, d4)) {
		if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();
		return 0;
	}

	const T bax = ((l3x*d2)-(l2x*d3));
	const T bay = ((l3y*d2)-(l2y*d3));
	const T baz = ((l3z*d2)-(l2z*d3));
	const T cax = ((l4x*d2)-(l2x*d4));
	const T cay = ((l4y*d2)-(l2y*d4));
	const T caz = ((l4z*d2)-(l2z*d4));
	const T qax = ((qx*d2)-l2x);
	const T qay = ((qy*d2)-l2y);
	const T qaz = ((qz*d2)-l2z);
	const T ba2 = ((ipow2(bax)+(ipow2(bay)+ipow2(baz)))*d4);
	const T ca2 = ((ipow2(cax)+(ipow2(cay)+ipow2(caz)))*d3);
	const T abcx = ((cax*ba2)-(bax*ca2));
	const T abcy = ((cay*ba2)-(bay*ca2));
	const T abcz = ((caz*ba2)-(baz*ca2));
	const T crossbcx = ((bay*caz)-(baz*cay));
	const T crossbcy = ((baz*cax)-(bax*caz));
	const T crossbcz = ((bax*cay)-(bay*cax));
	const T ccax = ((abcy*crossbcz)-(abcz*crossbcy));
	const T ccay = ((abcz*crossbcx)-(abcx*crossbcz));
	const T ccaz = ((abcx*crossbcy)-(abcy*crossbcx));
	const T c2 = ((ipow2(crossbcx)+(ipow2(crossbcy)+ipow2(crossbcz)))*(d3*d4));
	const T ret = ((qax*((qax*c2)-ccax))+((qay*((qay*c2)-ccay))+(qaz*((qaz*c2)-ccaz))));

	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();

	if constexpr (std::is_same<expansion, T>::value) {
		if (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW)) return INT_MAX;
	}
	return sgn(ret);
}

inline int inGabrielSphere_EIII(const genericPoint& a, const genericPoint& b, const genericPoint& c, const double& qx, const double& qy, const double& qz) {
	int ret;
	if ((ret = inGabrielSphere_EIII_t<interval_number, interval_number>(a, b, c, qx, qy, qz)) != 0) return ret;
	return inGabrielSphere_EIII_t<bigfloat, bigfloat>(a, b, c, qx, qy, qz);
}


template<class PT, class T> static inline int inGabrielSphere_IEEE_t(const genericPoint& q, const PT& ax, const PT& ay, const PT& az, const PT& bx, const PT& by, const PT& bz, const PT& cx, const PT& cy, const PT& cz) {
	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundUP();

	std::conditional_t<(std::is_same<expansion, T>::value), expansionPool, char> pool;

	if constexpr (std::is_same<expansion, T>::value) {
		expansion::initPool(&pool);
		feclearexcept(FE_ALL_EXCEPT);
	}

	T l1x, l1y, l1z, d1;
	if (!q.getLambda3D(l1x, l1y, l1z, d1)) {
		if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();
		return 0;
	}

	const T bax = (bx-ax);
	const T bay = (by-ay);
	const T baz = (bz-az);
	const T cax = (cx-ax);
	const T cay = (cy-ay);
	const T caz = (cz-az);
	const T crossbcx = ((bay*caz)-(baz*cay));
	const T crossbcy = ((baz*cax)-(bax*caz));
	const T crossbcz = ((bax*cay)-(bay*cax));
	const T ba2 = ((bax*bax)+((bay*bay)+(baz*baz)));
	const T ca2 = ((cax*cax)+((cay*cay)+(caz*caz)));
	const T calx = (cax*ba2);
	const T caly = (cay*ba2);
	const T calz = (caz*ba2);
	const T balx = (bax*ca2);
	const T baly = (bay*ca2);
	const T balz = (baz*ca2);
	const T abcx = (calx-balx);
	const T abcy = (caly-baly);
	const T abcz = (calz-balz);
	const T ccax = ((abcy*crossbcz)-(abcz*crossbcy));
	const T ccay = ((abcz*crossbcx)-(abcx*crossbcz));
	const T ccaz = ((abcx*crossbcy)-(abcy*crossbcx));
	const T c2 = (ipow2(crossbcx)+(ipow2(crossbcy)+ipow2(crossbcz)));
	const T qax = (l1x-(ax*d1));
	const T qay = (l1y-(ay*d1));
	const T qaz = (l1z-(az*d1));
	const T ret = ((qax*((qax*c2)-(ccax*d1)))+((qay*((qay*c2)-(ccay*d1)))+(qaz*((qaz*c2)-(ccaz*d1)))));

	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();

	if constexpr (std::is_same<expansion, T>::value) {
		if (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW)) return INT_MAX;
	}
	return sgn(ret);
}

inline int inGabrielSphere_IEEE(const genericPoint& q, const double& ax, const double& ay, const double& az, const double& bx, const double& by, const double& bz, const double& cx, const double& cy, const double& cz) {
	int ret;
	if ((ret = inGabrielSphere_IEEE_t<interval_number, interval_number>(q, ax, ay, az, bx, by, bz, cx, cy, cz)) != 0) return ret;
	return inGabrielSphere_IEEE_t<bigfloat, bigfloat>(q, ax, ay, az, bx, by, bz, cx, cy, cz);
}


template<class PT, class T> static inline int inGabrielSphere_IIEE_t(const genericPoint& q, const genericPoint& a, const PT& bx, const PT& by, const PT& bz, const PT& cx, const PT& cy, const PT& cz) {
	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundUP();

	std::conditional_t<(std::is_same<expansion, T>::value), expansionPool, char> pool;

	if constexpr (std::is_same<expansion, T>::value) {
		expansion::initPool(&pool);
		feclearexcept(FE_ALL_EXCEPT);
	}

	T l1x, l1y, l1z, d1;
	if (!q.getLambda3D(l1x, l1y, l1z, d1)) {
		if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();
		return 0;
	}

	T l2x, l2y, l2z, d2;
	if (!a.getLambda3D(l2x, l2y, l2z, d2)) {
		if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();
		return 0;
	}

	const T bax = ((bx*d2)-l2x);
	const T bay = ((by*d2)-l2y);
	const T baz = ((bz*d2)-l2z);
	const T cax = ((cx*d2)-l2x);
	const T cay = ((cy*d2)-l2y);
	const T caz = ((cz*d2)-l2z);
	const T crossbcx = ((bay*caz)-(baz*cay));
	const T crossbcy = ((baz*cax)-(bax*caz));
	const T crossbcz = ((bax*cay)-(bay*cax));
	const T ba2 = ((bax*bax)+((bay*bay)+(baz*baz)));
	const T ca2 = ((cax*cax)+((cay*cay)+(caz*caz)));
	const T calx = (cax*ba2);
	const T caly = (cay*ba2);
	const T calz = (caz*ba2);
	const T balx = (bax*ca2);
	const T baly = (bay*ca2);
	const T balz = (baz*ca2);
	const T abcx = (calx-balx);
	const T abcy = (caly-baly);
	const T abcz = (calz-balz);
	const T ccax = ((abcy*crossbcz)-(abcz*crossbcy));
	const T ccay = ((abcz*crossbcx)-(abcx*crossbcz));
	const T ccaz = ((abcx*crossbcy)-(abcy*crossbcx));
	const T c2 = (ipow2(crossbcx)+(ipow2(crossbcy)+ipow2(crossbcz)));
	const T qax = ((l1x*d2)-(l2x*d1));
	const T qay = ((l1y*d2)-(l2y*d1));
	const T qaz = ((l1z*d2)-(l2z*d1));
	const T ret = ((qax*((qax*c2)-(ccax*d1)))+((qay*((qay*c2)-(ccay*d1)))+(qaz*((qaz*c2)-(ccaz*d1)))));

	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();

	if constexpr (std::is_same<expansion, T>::value) {
		if (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW)) return INT_MAX;
	}
	return sgn(ret);
}

inline int inGabrielSphere_IIEE(const genericPoint& q, const genericPoint& a, const double& bx, const double& by, const double& bz, const double& cx, const double& cy, const double& cz) {
	int ret;
	if ((ret = inGabrielSphere_IIEE_t<interval_number, interval_number>(q, a, bx, by, bz, cx, cy, cz)) != 0) return ret;
	return inGabrielSphere_IIEE_t<bigfloat, bigfloat>(q, a, bx, by, bz, cx, cy, cz);
}


template<class PT, class T> static inline int inGabrielSphere_IIIE_t(const genericPoint& q, const genericPoint& a, const genericPoint& b, const PT& cx, const PT& cy, const PT& cz) {
	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundUP();

	std::conditional_t<(std::is_same<expansion, T>::value), expansionPool, char> pool;

	if constexpr (std::is_same<expansion, T>::value) {
		expansion::initPool(&pool);
		feclearexcept(FE_ALL_EXCEPT);
	}

	T l1x, l1y, l1z, d1;
	if (!q.getLambda3D(l1x, l1y, l1z, d1)) {
		if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();
		return 0;
	}

	T l2x, l2y, l2z, d2;
	if (!a.getLambda3D(l2x, l2y, l2z, d2)) {
		if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();
		return 0;
	}

	T l3x, l3y, l3z, d3;
	if (!b.getLambda3D(l3x, l3y, l3z, d3)) {
		if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();
		return 0;
	}

	const T bax = ((l3x*d2)-(l2x*d3));
	const T bay = ((l3y*d2)-(l2y*d3));
	const T baz = ((l3z*d2)-(l2z*d3));
	const T cax = ((cx*d2)-l2x);
	const T cay = ((cy*d2)-l2y);
	const T caz = ((cz*d2)-l2z);
	const T crossbcx = ((bay*caz)-(baz*cay));
	const T crossbcy = ((baz*cax)-(bax*caz));
	const T crossbcz = ((bax*cay)-(bay*cax));
	const T ba2 = ((bax*bax)+((bay*bay)+(baz*baz)));
	const T ca2 = (((cax*cax)+((cay*cay)+(caz*caz)))*d3);
	const T calx = (cax*ba2);
	const T caly = (cay*ba2);
	const T calz = (caz*ba2);
	const T balx = (bax*ca2);
	const T baly = (bay*ca2);
	const T balz = (baz*ca2);
	const T abcx = (calx-balx);
	const T abcy = (caly-baly);
	const T abcz = (calz-balz);
	const T ccax = ((abcy*crossbcz)-(abcz*crossbcy));
	const T ccay = ((abcz*crossbcx)-(abcx*crossbcz));
	const T ccaz = ((abcx*crossbcy)-(abcy*crossbcx));
	const T c2 = ((ipow2(crossbcx)+(ipow2(crossbcy)+ipow2(crossbcz)))*d3);
	const T qax = ((l1x*d2)-(l2x*d1));
	const T qay = ((l1y*d2)-(l2y*d1));
	const T qaz = ((l1z*d2)-(l2z*d1));
	const T ret = ((qax*((qax*c2)-(ccax*d1)))+((qay*((qay*c2)-(ccay*d1)))+(qaz*((qaz*c2)-(ccaz*d1)))));

	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();

	if constexpr (std::is_same<expansion, T>::value) {
		if (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW)) return INT_MAX;
	}
	return sgn(ret);
}

inline int inGabrielSphere_IIIE(const genericPoint& q, const genericPoint& a, const genericPoint& b, const double& cx, const double& cy, const double& cz) {
	int ret;
	if ((ret = inGabrielSphere_IIIE_t<interval_number, interval_number>(q, a, b, cx, cy, cz)) != 0) return ret;
	return inGabrielSphere_IIIE_t<bigfloat, bigfloat>(q, a, b, cx, cy, cz);
}


template<class PT, class T> static inline int inGabrielSphere_IIII_t(const genericPoint& q, const genericPoint& a, const genericPoint& b, const genericPoint& c) {
	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundUP();

	std::conditional_t<(std::is_same<expansion, T>::value), expansionPool, char> pool;

	if constexpr (std::is_same<expansion, T>::value) {
		expansion::initPool(&pool);
		feclearexcept(FE_ALL_EXCEPT);
	}

	T l1x, l1y, l1z, d1;
	if (!q.getLambda3D(l1x, l1y, l1z, d1)) {
		if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();
		return 0;
	}

	T l2x, l2y, l2z, d2;
	if (!a.getLambda3D(l2x, l2y, l2z, d2)) {
		if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();
		return 0;
	}

	T l3x, l3y, l3z, d3;
	if (!b.getLambda3D(l3x, l3y, l3z, d3)) {
		if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();
		return 0;
	}

	T l4x, l4y, l4z, d4;
	if (!c.getLambda3D(l4x, l4y, l4z, d4)) {
		if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();
		return 0;
	}

	const T bax = ((l3x*d2)-(l2x*d3));
	const T bay = ((l3y*d2)-(l2y*d3));
	const T baz = ((l3z*d2)-(l2z*d3));
	const T cax = ((l4x*d2)-(l2x*d4));
	const T cay = ((l4y*d2)-(l2y*d4));
	const T caz = ((l4z*d2)-(l2z*d4));
	const T crossbcx = ((bay*caz)-(baz*cay));
	const T crossbcy = ((baz*cax)-(bax*caz));
	const T crossbcz = ((bax*cay)-(bay*cax));
	const T ba2 = (((bax*bax)+((bay*bay)+(baz*baz)))*d4);
	const T ca2 = (((cax*cax)+((cay*cay)+(caz*caz)))*d3);
	const T calx = (cax*ba2);
	const T caly = (cay*ba2);
	const T calz = (caz*ba2);
	const T balx = (bax*ca2);
	const T baly = (bay*ca2);
	const T balz = (baz*ca2);
	const T abcx = (calx-balx);
	const T abcy = (caly-baly);
	const T abcz = (calz-balz);
	const T ccax = ((abcy*crossbcz)-(abcz*crossbcy));
	const T ccay = ((abcz*crossbcx)-(abcx*crossbcz));
	const T ccaz = ((abcx*crossbcy)-(abcy*crossbcx));
	const T c2 = ((ipow2(crossbcx)+(ipow2(crossbcy)+ipow2(crossbcz)))*(d3*d4));
	const T qax = ((l1x*d2)-(l2x*d1));
	const T qay = ((l1y*d2)-(l2y*d1));
	const T qaz = ((l1z*d2)-(l2z*d1));
	const T ret = ((qax*((qax*c2)-(ccax*d1)))+((qay*((qay*c2)-(ccay*d1)))+(qaz*((qaz*c2)-(ccaz*d1)))));

	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();

	if constexpr (std::is_same<expansion, T>::value) {
		if (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW)) return INT_MAX;
	}
	return sgn(ret);
}

inline int inGabrielSphere_IIII(const genericPoint& q, const genericPoint& a, const genericPoint& b, const genericPoint& c) {
	int ret;
	if ((ret = inGabrielSphere_IIII_t<interval_number, interval_number>(q, a, b, c)) != 0) return ret;
	return inGabrielSphere_IIII_t<bigfloat, bigfloat>(q, a, b, c);
}


template<class PT, class T> static inline int inSphere_IEEEE_t(const genericPoint& p1, const PT& pbx, const PT& pby, const PT& pbz, const PT& pcx, const PT& pcy, const PT& pcz, const PT& pdx, const PT& pdy, const PT& pdz, const PT& pex, const PT& pey, const PT& pez) {
	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundUP();

	std::conditional_t<(std::is_same<expansion, T>::value), expansionPool, char> pool;

	if constexpr (std::is_same<expansion, T>::value) {
		expansion::initPool(&pool);
		feclearexcept(FE_ALL_EXCEPT);
	}

	T l1x, l1y, l1z, d1;
	if (!p1.getLambda3D(l1x, l1y, l1z, d1)) {
		if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();
		return 0;
	}

	const T aex = (l1x-(pex*d1));
	const T aey = (l1y-(pey*d1));
	const T aez = (l1z-(pez*d1));
	const T alift = (ipow2(aex)+(ipow2(aey)+ipow2(aez)));
	const T bex = (pbx-pex);
	const T bey = (pby-pey);
	const T bez = (pbz-pez);
	const T blift = (ipow2(bex)+(ipow2(bey)+ipow2(bez)));
	const T cex = (pcx-pex);
	const T cey = (pcy-pey);
	const T cez = (pcz-pez);
	const T clift = (ipow2(cex)+(ipow2(cey)+ipow2(cez)));
	const T dex = (pdx-pex);
	const T dey = (pdy-pey);
	const T dez = (pdz-pez);
	const T dlift = (ipow2(dex)+(ipow2(dey)+ipow2(dez)));
	const T ab = ((aex*bey)-(bex*aey));
	const T bc = ((bex*cey)-(cex*bey));
	const T cd = ((cex*dey)-(dex*cey));
	const T da = ((dex*aey)-(aex*dey));
	const T ac = ((aex*cey)-(cex*aey));
	const T bd = ((bex*dey)-(dex*bey));
	const T abc = (((aez*bc)-(bez*ac))+(cez*ab));
	const T bcd = (((bez*cd)-(cez*bd))+(dez*bc));
	const T cda = ((cez*da)+((dez*ac)+(aez*cd)));
	const T dab = ((dez*ab)+((aez*bd)+(bez*da)));
	const T d = (((clift*dab)-(dlift*abc))+((alift*bcd)-(blift*cda)));

	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();

	if constexpr (std::is_same<expansion, T>::value) {
		if (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW)) return INT_MAX;
	}
	return sgn(d);
}

inline int inSphere_IEEEE(const genericPoint& p1, const double& pbx, const double& pby, const double& pbz, const double& pcx, const double& pcy, const double& pcz, const double& pdx, const double& pdy, const double& pdz, const double& pex, const double& pey, const double& pez) {
	int ret;
	if ((ret = inSphere_IEEEE_t<interval_number, interval_number>(p1, pbx, pby, pbz, pcx, pcy, pcz, pdx, pdy, pdz, pex, pey, pez)) != 0) return ret;
	return inSphere_IEEEE_t<bigfloat, bigfloat>(p1, pbx, pby, pbz, pcx, pcy, pcz, pdx, pdy, pdz, pex, pey, pez);
}


template<class PT, class T> static inline int inSphere_IIEEE_t(const genericPoint& p1, const genericPoint& p2, const PT& pcx, const PT& pcy, const PT& pcz, const PT& pdx, const PT& pdy, const PT& pdz, const PT& pex, const PT& pey, const PT& pez) {
	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundUP();

	std::conditional_t<(std::is_same<expansion, T>::value), expansionPool, char> pool;

	if constexpr (std::is_same<expansion, T>::value) {
		expansion::initPool(&pool);
		feclearexcept(FE_ALL_EXCEPT);
	}

	T l1x, l1y, l1z, d1;
	if (!p1.getLambda3D(l1x, l1y, l1z, d1)) {
		if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();
		return 0;
	}

	T l2x, l2y, l2z, d2;
	if (!p2.getLambda3D(l2x, l2y, l2z, d2)) {
		if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();
		return 0;
	}

	const T aex = (l1x-(pex*d1));
	const T aey = (l1y-(pey*d1));
	const T aez = (l1z-(pez*d1));
	const T alift = (ipow2(aex)+(ipow2(aey)+ipow2(aez)));
	const T bex = (l2x-(pex*d2));
	const T bey = (l2y-(pey*d2));
	const T bez = (l2z-(pez*d2));
	const T blift = (ipow2(bex)+(ipow2(bey)+ipow2(bez)));
	const T cex = (pcx-pex);
	const T cey = (pcy-pey);
	const T cez = (pcz-pez);
	const T clift = (ipow2(cex)+(ipow2(cey)+ipow2(cez)));
	const T dex = (pdx-pex);
	const T dey = (pdy-pey);
	const T dez = (pdz-pez);
	const T dlift = (ipow2(dex)+(ipow2(dey)+ipow2(dez)));
	const T ab = ((aex*bey)-(bex*aey));
	const T bc = ((bex*cey)-(cex*bey));
	const T cd = ((cex*dey)-(dex*cey));
	const T da = ((dex*aey)-(aex*dey));
	const T ac = ((aex*cey)-(cex*aey));
	const T bd = ((bex*dey)-(dex*bey));
	const T abc = (((aez*bc)-(bez*ac))+(cez*ab));
	const T bcd = (((bez*cd)-(cez*bd))+(dez*bc));
	const T cda = ((cez*da)+((dez*ac)+(aez*cd)));
	const T dab = ((dez*ab)+((aez*bd)+(bez*da)));
	const T d = (((clift*dab)-(dlift*abc))+((alift*bcd)-(blift*cda)));

	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();

	if constexpr (std::is_same<expansion, T>::value) {
		if (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW)) return INT_MAX;
	}
	return sgn(d);
}

inline int inSphere_IIEEE(const genericPoint& p1, const genericPoint& p2, const double& pcx, const double& pcy, const double& pcz, const double& pdx, const double& pdy, const double& pdz, const double& pex, const double& pey, const double& pez) {
	int ret;
	if ((ret = inSphere_IIEEE_t<interval_number, interval_number>(p1, p2, pcx, pcy, pcz, pdx, pdy, pdz, pex, pey, pez)) != 0) return ret;
	return inSphere_IIEEE_t<bigfloat, bigfloat>(p1, p2, pcx, pcy, pcz, pdx, pdy, pdz, pex, pey, pez);
}


template<class PT, class T> static inline int inSphere_IIIEE_t(const genericPoint& p1, const genericPoint& p2, const genericPoint& p3, const PT& pdx, const PT& pdy, const PT& pdz, const PT& pex, const PT& pey, const PT& pez) {
	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundUP();

	std::conditional_t<(std::is_same<expansion, T>::value), expansionPool, char> pool;

	if constexpr (std::is_same<expansion, T>::value) {
		expansion::initPool(&pool);
		feclearexcept(FE_ALL_EXCEPT);
	}

	T l1x, l1y, l1z, d1;
	if (!p1.getLambda3D(l1x, l1y, l1z, d1)) {
		if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();
		return 0;
	}

	T l2x, l2y, l2z, d2;
	if (!p2.getLambda3D(l2x, l2y, l2z, d2)) {
		if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();
		return 0;
	}

	T l3x, l3y, l3z, d3;
	if (!p3.getLambda3D(l3x, l3y, l3z, d3)) {
		if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();
		return 0;
	}

	const T aex = (l1x-(pex*d1));
	const T aey = (l1y-(pey*d1));
	const T aez = (l1z-(pez*d1));
	const T alift = (ipow2(aex)+(ipow2(aey)+ipow2(aez)));
	const T bex = (l2x-(pex*d2));
	const T bey = (l2y-(pey*d2));
	const T bez = (l2z-(pez*d2));
	const T blift = (ipow2(bex)+(ipow2(bey)+ipow2(bez)));
	const T cex = (l3x-(pex*d3));
	const T cey = (l3y-(pey*d3));
	const T cez = (l3z-(pez*d3));
	const T clift = (ipow2(cex)+(ipow2(cey)+ipow2(cez)));
	const T dex = (pdx-pex);
	const T dey = (pdy-pey);
	const T dez = (pdz-pez);
	const T dlift = (ipow2(dex)+(ipow2(dey)+ipow2(dez)));
	const T ab = ((aex*bey)-(bex*aey));
	const T bc = ((bex*cey)-(cex*bey));
	const T cd = ((cex*dey)-(dex*cey));
	const T da = ((dex*aey)-(aex*dey));
	const T ac = ((aex*cey)-(cex*aey));
	const T bd = ((bex*dey)-(dex*bey));
	const T abc = (((aez*bc)-(bez*ac))+(cez*ab));
	const T bcd = (((bez*cd)-(cez*bd))+(dez*bc));
	const T cda = ((cez*da)+((dez*ac)+(aez*cd)));
	const T dab = ((dez*ab)+((aez*bd)+(bez*da)));
	const T d = (((clift*dab)-(dlift*abc))+((alift*bcd)-(blift*cda)));

	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();

	if constexpr (std::is_same<expansion, T>::value) {
		if (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW)) return INT_MAX;
	}
	return sgn(d);
}

inline int inSphere_IIIEE(const genericPoint& p1, const genericPoint& p2, const genericPoint& p3, const double& pdx, const double& pdy, const double& pdz, const double& pex, const double& pey, const double& pez) {
	int ret;
	if ((ret = inSphere_IIIEE_t<interval_number, interval_number>(p1, p2, p3, pdx, pdy, pdz, pex, pey, pez)) != 0) return ret;
	return inSphere_IIIEE_t<bigfloat, bigfloat>(p1, p2, p3, pdx, pdy, pdz, pex, pey, pez);
}


template<class PT, class T> static inline int inSphere_IIIIE_t(const genericPoint& p1, const genericPoint& p2, const genericPoint& p3, const genericPoint& p4, const PT& pex, const PT& pey, const PT& pez) {
	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundUP();

	std::conditional_t<(std::is_same<expansion, T>::value), expansionPool, char> pool;

	if constexpr (std::is_same<expansion, T>::value) {
		expansion::initPool(&pool);
		feclearexcept(FE_ALL_EXCEPT);
	}

	T l1x, l1y, l1z, d1;
	if (!p1.getLambda3D(l1x, l1y, l1z, d1)) {
		if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();
		return 0;
	}

	T l2x, l2y, l2z, d2;
	if (!p2.getLambda3D(l2x, l2y, l2z, d2)) {
		if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();
		return 0;
	}

	T l3x, l3y, l3z, d3;
	if (!p3.getLambda3D(l3x, l3y, l3z, d3)) {
		if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();
		return 0;
	}

	T l4x, l4y, l4z, d4;
	if (!p4.getLambda3D(l4x, l4y, l4z, d4)) {
		if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();
		return 0;
	}

	const T aex = (l1x-(pex*d1));
	const T aey = (l1y-(pey*d1));
	const T aez = (l1z-(pez*d1));
	const T alift = (ipow2(aex)+(ipow2(aey)+ipow2(aez)));
	const T bex = (l2x-(pex*d2));
	const T bey = (l2y-(pey*d2));
	const T bez = (l2z-(pez*d2));
	const T blift = (ipow2(bex)+(ipow2(bey)+ipow2(bez)));
	const T cex = (l3x-(pex*d3));
	const T cey = (l3y-(pey*d3));
	const T cez = (l3z-(pez*d3));
	const T clift = (ipow2(cex)+(ipow2(cey)+ipow2(cez)));
	const T dex = (l4x-(pex*d4));
	const T dey = (l4y-(pey*d4));
	const T dez = (l4z-(pez*d4));
	const T dlift = (ipow2(dex)+(ipow2(dey)+ipow2(dez)));
	const T ab = ((aex*bey)-(bex*aey));
	const T bc = ((bex*cey)-(cex*bey));
	const T cd = ((cex*dey)-(dex*cey));
	const T da = ((dex*aey)-(aex*dey));
	const T ac = ((aex*cey)-(cex*aey));
	const T bd = ((bex*dey)-(dex*bey));
	const T abc = (((aez*bc)-(bez*ac))+(cez*ab));
	const T bcd = (((bez*cd)-(cez*bd))+(dez*bc));
	const T cda = ((cez*da)+((dez*ac)+(aez*cd)));
	const T dab = ((dez*ab)+((aez*bd)+(bez*da)));
	const T d = (((clift*dab)-(dlift*abc))+((alift*bcd)-(blift*cda)));

	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();

	if constexpr (std::is_same<expansion, T>::value) {
		if (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW)) return INT_MAX;
	}
	return sgn(d);
}

inline int inSphere_IIIIE(const genericPoint& p1, const genericPoint& p2, const genericPoint& p3, const genericPoint& p4, const double& pex, const double& pey, const double& pez) {
	int ret;
	if ((ret = inSphere_IIIIE_t<interval_number, interval_number>(p1, p2, p3, p4, pex, pey, pez)) != 0) return ret;
	return inSphere_IIIIE_t<bigfloat, bigfloat>(p1, p2, p3, p4, pex, pey, pez);
}


template<class PT, class T> static inline int inSphere_IIIII_t(const genericPoint& p1, const genericPoint& p2, const genericPoint& p3, const genericPoint& p4, const genericPoint& p5) {
	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundUP();

	std::conditional_t<(std::is_same<expansion, T>::value), expansionPool, char> pool;

	if constexpr (std::is_same<expansion, T>::value) {
		expansion::initPool(&pool);
		feclearexcept(FE_ALL_EXCEPT);
	}

	T l1x, l1y, l1z, d1;
	if (!p1.getLambda3D(l1x, l1y, l1z, d1)) {
		if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();
		return 0;
	}

	T l2x, l2y, l2z, d2;
	if (!p2.getLambda3D(l2x, l2y, l2z, d2)) {
		if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();
		return 0;
	}

	T l3x, l3y, l3z, d3;
	if (!p3.getLambda3D(l3x, l3y, l3z, d3)) {
		if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();
		return 0;
	}

	T l4x, l4y, l4z, d4;
	if (!p4.getLambda3D(l4x, l4y, l4z, d4)) {
		if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();
		return 0;
	}

	T l5x, l5y, l5z, d5;
	if (!p5.getLambda3D(l5x, l5y, l5z, d5)) {
		if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();
		return 0;
	}

	const T aex = ((l1x*d5)-(l5x*d1));
	const T aey = ((l1y*d5)-(l5y*d1));
	const T aez = ((l1z*d5)-(l5z*d1));
	const T alift = (ipow2(aex)+(ipow2(aey)+ipow2(aez)));
	const T bex = ((l2x*d5)-(l5x*d2));
	const T bey = ((l2y*d5)-(l5y*d2));
	const T bez = ((l2z*d5)-(l5z*d2));
	const T blift = (ipow2(bex)+(ipow2(bey)+ipow2(bez)));
	const T cex = ((l3x*d5)-(l5x*d3));
	const T cey = ((l3y*d5)-(l5y*d3));
	const T cez = ((l3z*d5)-(l5z*d3));
	const T clift = (ipow2(cex)+(ipow2(cey)+ipow2(cez)));
	const T dex = ((l4x*d5)-(l5x*d4));
	const T dey = ((l4y*d5)-(l5y*d4));
	const T dez = ((l4z*d5)-(l5z*d4));
	const T dlift = (ipow2(dex)+(ipow2(dey)+ipow2(dez)));
	const T ab = ((aex*bey)-(bex*aey));
	const T bc = ((bex*cey)-(cex*bey));
	const T cd = ((cex*dey)-(dex*cey));
	const T da = ((dex*aey)-(aex*dey));
	const T ac = ((aex*cey)-(cex*aey));
	const T bd = ((bex*dey)-(dex*bey));
	const T abc = (((aez*bc)-(bez*ac))+(cez*ab));
	const T bcd = (((bez*cd)-(cez*bd))+(dez*bc));
	const T cda = ((cez*da)+((dez*ac)+(aez*cd)));
	const T dab = ((dez*ab)+((aez*bd)+(bez*da)));
	const T d = (((clift*dab)-(dlift*abc))+((alift*bcd)-(blift*cda)));

	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();

	if constexpr (std::is_same<expansion, T>::value) {
		if (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW)) return INT_MAX;
	}
	return sgn(d);
}

inline int inSphere_IIIII(const genericPoint& p1, const genericPoint& p2, const genericPoint& p3, const genericPoint& p4, const genericPoint& p5) {
	int ret;
	if ((ret = inSphere_IIIII_t<interval_number, interval_number>(p1, p2, p3, p4, p5)) != 0) return ret;
	return inSphere_IIIII_t<bigfloat, bigfloat>(p1, p2, p3, p4, p5);
}


template<class PT, class T> static inline int lessThanOnX_IE_t(const genericPoint& p1, const PT& bx) {
	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundUP();

	std::conditional_t<(std::is_same<expansion, T>::value), expansionPool, char> pool;

	if constexpr (std::is_same<expansion, T>::value) {
		expansion::initPool(&pool);
		feclearexcept(FE_ALL_EXCEPT);
	}

	T l1x, l1y, l1z, d1;
	if (!p1.getLambda3D(l1x, l1y, l1z, d1)) {
		if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();
		return 0;
	}

	const T dbx = (bx*d1);
	const T kx = (l1x-dbx);

	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();

	if constexpr (std::is_same<expansion, T>::value) {
		if (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW)) return INT_MAX;
	}
	return sgn(kx);
}

inline int lessThanOnX_IE(const genericPoint& p1, const double& bx) {
	int ret;
	if ((ret = lessThanOnX_IE_t<interval_number, interval_number>(p1, bx)) != 0) return ret;
	return lessThanOnX_IE_t<bigfloat, bigfloat>(p1, bx);
}


template<class PT, class T> static inline int lessThanOnX_II_t(const genericPoint& p1, const genericPoint& p2) {
	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundUP();

	std::conditional_t<(std::is_same<expansion, T>::value), expansionPool, char> pool;

	if constexpr (std::is_same<expansion, T>::value) {
		expansion::initPool(&pool);
		feclearexcept(FE_ALL_EXCEPT);
	}

	T l1x, l1y, l1z, d1;
	if (!p1.getLambda3D(l1x, l1y, l1z, d1)) {
		if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();
		return 0;
	}

	T l2x, l2y, l2z, d2;
	if (!p2.getLambda3D(l2x, l2y, l2z, d2)) {
		if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();
		return 0;
	}

	const T k1 = (d2*l1x);
	const T k2 = (d1*l2x);
	const T kx = (k1-k2);

	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();

	if constexpr (std::is_same<expansion, T>::value) {
		if (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW)) return INT_MAX;
	}
	return sgn(kx);
}

inline int lessThanOnX_II(const genericPoint& p1, const genericPoint& p2) {
	int ret;
	if ((ret = lessThanOnX_II_t<interval_number, interval_number>(p1, p2)) != 0) return ret;
	return lessThanOnX_II_t<bigfloat, bigfloat>(p1, p2);
}


template<class PT, class T> static inline int lessThanOnY_IE_t(const genericPoint& p1, const PT& by) {
	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundUP();

	std::conditional_t<(std::is_same<expansion, T>::value), expansionPool, char> pool;

	if constexpr (std::is_same<expansion, T>::value) {
		expansion::initPool(&pool);
		feclearexcept(FE_ALL_EXCEPT);
	}

	T l1x, l1y, l1z, d1;
	if (!p1.getLambda3D(l1x, l1y, l1z, d1)) {
		if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();
		return 0;
	}

	const T dby = (by*d1);
	const T ky = (l1y-dby);

	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();

	if constexpr (std::is_same<expansion, T>::value) {
		if (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW)) return INT_MAX;
	}
	return sgn(ky);
}

inline int lessThanOnY_IE(const genericPoint& p1, const double& by) {
	int ret;
	if ((ret = lessThanOnY_IE_t<interval_number, interval_number>(p1, by)) != 0) return ret;
	return lessThanOnY_IE_t<bigfloat, bigfloat>(p1, by);
}


template<class PT, class T> static inline int lessThanOnY_II_t(const genericPoint& p1, const genericPoint& p2) {
	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundUP();

	std::conditional_t<(std::is_same<expansion, T>::value), expansionPool, char> pool;

	if constexpr (std::is_same<expansion, T>::value) {
		expansion::initPool(&pool);
		feclearexcept(FE_ALL_EXCEPT);
	}

	T l1x, l1y, l1z, d1;
	if (!p1.getLambda3D(l1x, l1y, l1z, d1)) {
		if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();
		return 0;
	}

	T l2x, l2y, l2z, d2;
	if (!p2.getLambda3D(l2x, l2y, l2z, d2)) {
		if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();
		return 0;
	}

	const T k1 = (d2*l1y);
	const T k2 = (d1*l2y);
	const T ky = (k1-k2);

	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();

	if constexpr (std::is_same<expansion, T>::value) {
		if (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW)) return INT_MAX;
	}
	return sgn(ky);
}

inline int lessThanOnY_II(const genericPoint& p1, const genericPoint& p2) {
	int ret;
	if ((ret = lessThanOnY_II_t<interval_number, interval_number>(p1, p2)) != 0) return ret;
	return lessThanOnY_II_t<bigfloat, bigfloat>(p1, p2);
}


template<class PT, class T> static inline int lessThanOnZ_IE_t(const genericPoint& p1, const PT& bz) {
	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundUP();

	std::conditional_t<(std::is_same<expansion, T>::value), expansionPool, char> pool;

	if constexpr (std::is_same<expansion, T>::value) {
		expansion::initPool(&pool);
		feclearexcept(FE_ALL_EXCEPT);
	}

	T l1x, l1y, l1z, d1;
	if (!p1.getLambda3D(l1x, l1y, l1z, d1)) {
		if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();
		return 0;
	}

	const T dbz = (bz*d1);
	const T kz = (l1z-dbz);

	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();

	if constexpr (std::is_same<expansion, T>::value) {
		if (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW)) return INT_MAX;
	}
	return sgn(kz);
}

inline int lessThanOnZ_IE(const genericPoint& p1, const double& bz) {
	int ret;
	if ((ret = lessThanOnZ_IE_t<interval_number, interval_number>(p1, bz)) != 0) return ret;
	return lessThanOnZ_IE_t<bigfloat, bigfloat>(p1, bz);
}


template<class PT, class T> static inline int lessThanOnZ_II_t(const genericPoint& p1, const genericPoint& p2) {
	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundUP();

	std::conditional_t<(std::is_same<expansion, T>::value), expansionPool, char> pool;

	if constexpr (std::is_same<expansion, T>::value) {
		expansion::initPool(&pool);
		feclearexcept(FE_ALL_EXCEPT);
	}

	T l1x, l1y, l1z, d1;
	if (!p1.getLambda3D(l1x, l1y, l1z, d1)) {
		if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();
		return 0;
	}

	T l2x, l2y, l2z, d2;
	if (!p2.getLambda3D(l2x, l2y, l2z, d2)) {
		if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();
		return 0;
	}

	const T k1 = (d2*l1z);
	const T k2 = (d1*l2z);
	const T kz = (k1-k2);

	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();

	if constexpr (std::is_same<expansion, T>::value) {
		if (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW)) return INT_MAX;
	}
	return sgn(kz);
}

inline int lessThanOnZ_II(const genericPoint& p1, const genericPoint& p2) {
	int ret;
	if ((ret = lessThanOnZ_II_t<interval_number, interval_number>(p1, p2)) != 0) return ret;
	return lessThanOnZ_II_t<bigfloat, bigfloat>(p1, p2);
}


template<class PT, class T> static inline int orient2dxy_indirect_IEE_t(const genericPoint& p1, const PT& p2x, const PT& p2y, const PT& p3x, const PT& p3y) {
	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundUP();

	std::conditional_t<(std::is_same<expansion, T>::value), expansionPool, char> pool;

	if constexpr (std::is_same<expansion, T>::value) {
		expansion::initPool(&pool);
		feclearexcept(FE_ALL_EXCEPT);
	}

	T l1x, l1y, l1z, d1;
	if (!p1.getLambda3D(l1x, l1y, l1z, d1)) {
		if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();
		return 0;
	}

	const T t1x = (p2y-p3y);
	const T t1y = (p3x-p2x);
	const T e2 = (l1x*t1x);
	const T e3 = (l1y*t1y);
	const T e = (e2+e3);
	const T pr1 = (p2x*p3y);
	const T pr2 = (p2y*p3x);
	const T pr = (pr1-pr2);
	const T dpr = (d1*pr);
	const T det = (dpr+e);

	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();

	if constexpr (std::is_same<expansion, T>::value) {
		if (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW)) return INT_MAX;
	}
	return sgn(det);
}

inline int orient2dxy_indirect_IEE(const genericPoint& p1, const double& p2x, const double& p2y, const double& p3x, const double& p3y) {
	int ret;
	if ((ret = orient2dxy_indirect_IEE_t<interval_number, interval_number>(p1, p2x, p2y, p3x, p3y)) != 0) return ret;
	return orient2dxy_indirect_IEE_t<bigfloat, bigfloat>(p1, p2x, p2y, p3x, p3y);
}


template<class PT, class T> static inline int orient2dxy_indirect_IIE_t(const genericPoint& p1, const genericPoint& p2, const PT& op3x, const PT& op3y) {
	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundUP();

	std::conditional_t<(std::is_same<expansion, T>::value), expansionPool, char> pool;

	if constexpr (std::is_same<expansion, T>::value) {
		expansion::initPool(&pool);
		feclearexcept(FE_ALL_EXCEPT);
	}

	T l1x, l1y, l1z, d1;
	if (!p1.getLambda3D(l1x, l1y, l1z, d1)) {
		if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();
		return 0;
	}

	T l2x, l2y, l2z, d2;
	if (!p2.getLambda3D(l2x, l2y, l2z, d2)) {
		if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();
		return 0;
	}

	const T a = (d1*l2x);
	const T b = (d2*l1x);
	const T c = (d1*op3y);
	const T e = (d1*l2y);
	const T f = (d2*l1y);
	const T g = (d1*op3x);
	const T ab = (a-b);
	const T cd = (c-l1y);
	const T ef = (e-f);
	const T gh = (g-l1x);
	const T abcd = (ab*cd);
	const T efgh = (ef*gh);
	const T L = (abcd-efgh);

	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();

	if constexpr (std::is_same<expansion, T>::value) {
		if (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW)) return INT_MAX;
	}
	return sgn(L);
}

inline int orient2dxy_indirect_IIE(const genericPoint& p1, const genericPoint& p2, const double& op3x, const double& op3y) {
	int ret;
	if ((ret = orient2dxy_indirect_IIE_t<interval_number, interval_number>(p1, p2, op3x, op3y)) != 0) return ret;
	return orient2dxy_indirect_IIE_t<bigfloat, bigfloat>(p1, p2, op3x, op3y);
}


template<class PT, class T> static inline int orient2dxy_indirect_III_t(const genericPoint& p1, const genericPoint& p2, const genericPoint& p3) {
	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundUP();

	std::conditional_t<(std::is_same<expansion, T>::value), expansionPool, char> pool;

	if constexpr (std::is_same<expansion, T>::value) {
		expansion::initPool(&pool);
		feclearexcept(FE_ALL_EXCEPT);
	}

	T l1x, l1y, l1z, d1;
	if (!p1.getLambda3D(l1x, l1y, l1z, d1)) {
		if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();
		return 0;
	}

	T l2x, l2y, l2z, d2;
	if (!p2.getLambda3D(l2x, l2y, l2z, d2)) {
		if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();
		return 0;
	}

	T l3x, l3y, l3z, d3;
	if (!p3.getLambda3D(l3x, l3y, l3z, d3)) {
		if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();
		return 0;
	}

	const T a = (d1*l2x);
	const T b = (d2*l1x);
	const T c = (d1*l3y);
	const T d = (d3*l1y);
	const T e = (d1*l2y);
	const T f = (d2*l1y);
	const T g = (d1*l3x);
	const T h = (d3*l1x);
	const T ab = (a-b);
	const T cd = (c-d);
	const T ef = (e-f);
	const T gh = (g-h);
	const T abcd = (ab*cd);
	const T efgh = (ef*gh);
	const T L = (abcd-efgh);

	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();

	if constexpr (std::is_same<expansion, T>::value) {
		if (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW)) return INT_MAX;
	}
	return sgn(L);
}

inline int orient2dxy_indirect_III(const genericPoint& p1, const genericPoint& p2, const genericPoint& p3) {
	int ret;
	if ((ret = orient2dxy_indirect_III_t<interval_number, interval_number>(p1, p2, p3)) != 0) return ret;
	return orient2dxy_indirect_III_t<bigfloat, bigfloat>(p1, p2, p3);
}


template<class PT, class T> static inline int orient2dyz_indirect_IEE_t(const genericPoint& p1, const PT& p2x, const PT& p2y, const PT& p3x, const PT& p3y) {
	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundUP();

	std::conditional_t<(std::is_same<expansion, T>::value), expansionPool, char> pool;

	if constexpr (std::is_same<expansion, T>::value) {
		expansion::initPool(&pool);
		feclearexcept(FE_ALL_EXCEPT);
	}

	T l1z, l1x, l1y, d1;
	if (!p1.getLambda3D(l1z, l1x, l1y, d1)) {
		if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();
		return 0;
	}

	const T t1x = (p2y-p3y);
	const T t1y = (p3x-p2x);
	const T e2 = (l1x*t1x);
	const T e3 = (l1y*t1y);
	const T e = (e2+e3);
	const T pr1 = (p2x*p3y);
	const T pr2 = (p2y*p3x);
	const T pr = (pr1-pr2);
	const T dpr = (d1*pr);
	const T det = (dpr+e);

	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();

	if constexpr (std::is_same<expansion, T>::value) {
		if (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW)) return INT_MAX;
	}
	return sgn(det);
}

inline int orient2dyz_indirect_IEE(const genericPoint& p1, const double& p2x, const double& p2y, const double& p3x, const double& p3y) {
	int ret;
	if ((ret = orient2dyz_indirect_IEE_t<interval_number, interval_number>(p1, p2x, p2y, p3x, p3y)) != 0) return ret;
	return orient2dyz_indirect_IEE_t<bigfloat, bigfloat>(p1, p2x, p2y, p3x, p3y);
}


template<class PT, class T> static inline int orient2dyz_indirect_IIE_t(const genericPoint& p1, const genericPoint& p2, const PT& op3x, const PT& op3y) {
	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundUP();

	std::conditional_t<(std::is_same<expansion, T>::value), expansionPool, char> pool;

	if constexpr (std::is_same<expansion, T>::value) {
		expansion::initPool(&pool);
		feclearexcept(FE_ALL_EXCEPT);
	}

	T l1z, l1x, l1y, d1;
	if (!p1.getLambda3D(l1z, l1x, l1y, d1)) {
		if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();
		return 0;
	}

	T l2z, l2x, l2y, d2;
	if (!p2.getLambda3D(l2z, l2x, l2y, d2)) {
		if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();
		return 0;
	}

	const T a = (d1*l2x);
	const T b = (d2*l1x);
	const T c = (d1*op3y);
	const T e = (d1*l2y);
	const T f = (d2*l1y);
	const T g = (d1*op3x);
	const T ab = (a-b);
	const T cd = (c-l1y);
	const T ef = (e-f);
	const T gh = (g-l1x);
	const T abcd = (ab*cd);
	const T efgh = (ef*gh);
	const T L = (abcd-efgh);

	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();

	if constexpr (std::is_same<expansion, T>::value) {
		if (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW)) return INT_MAX;
	}
	return sgn(L);
}

inline int orient2dyz_indirect_IIE(const genericPoint& p1, const genericPoint& p2, const double& op3x, const double& op3y) {
	int ret;
	if ((ret = orient2dyz_indirect_IIE_t<interval_number, interval_number>(p1, p2, op3x, op3y)) != 0) return ret;
	return orient2dyz_indirect_IIE_t<bigfloat, bigfloat>(p1, p2, op3x, op3y);
}


template<class PT, class T> static inline int orient2dyz_indirect_III_t(const genericPoint& p1, const genericPoint& p2, const genericPoint& p3) {
	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundUP();

	std::conditional_t<(std::is_same<expansion, T>::value), expansionPool, char> pool;

	if constexpr (std::is_same<expansion, T>::value) {
		expansion::initPool(&pool);
		feclearexcept(FE_ALL_EXCEPT);
	}

	T l1z, l1x, l1y, d1;
	if (!p1.getLambda3D(l1z, l1x, l1y, d1)) {
		if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();
		return 0;
	}

	T l2z, l2x, l2y, d2;
	if (!p2.getLambda3D(l2z, l2x, l2y, d2)) {
		if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();
		return 0;
	}

	T l3z, l3x, l3y, d3;
	if (!p3.getLambda3D(l3z, l3x, l3y, d3)) {
		if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();
		return 0;
	}

	const T a = (d1*l2x);
	const T b = (d2*l1x);
	const T c = (d1*l3y);
	const T d = (d3*l1y);
	const T e = (d1*l2y);
	const T f = (d2*l1y);
	const T g = (d1*l3x);
	const T h = (d3*l1x);
	const T ab = (a-b);
	const T cd = (c-d);
	const T ef = (e-f);
	const T gh = (g-h);
	const T abcd = (ab*cd);
	const T efgh = (ef*gh);
	const T L = (abcd-efgh);

	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();

	if constexpr (std::is_same<expansion, T>::value) {
		if (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW)) return INT_MAX;
	}
	return sgn(L);
}

inline int orient2dyz_indirect_III(const genericPoint& p1, const genericPoint& p2, const genericPoint& p3) {
	int ret;
	if ((ret = orient2dyz_indirect_III_t<interval_number, interval_number>(p1, p2, p3)) != 0) return ret;
	return orient2dyz_indirect_III_t<bigfloat, bigfloat>(p1, p2, p3);
}


template<class PT, class T> static inline int orient2dzx_indirect_IEE_t(const genericPoint& p1, const PT& p2x, const PT& p2y, const PT& p3x, const PT& p3y) {
	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundUP();

	std::conditional_t<(std::is_same<expansion, T>::value), expansionPool, char> pool;

	if constexpr (std::is_same<expansion, T>::value) {
		expansion::initPool(&pool);
		feclearexcept(FE_ALL_EXCEPT);
	}

	T l1y, l1z, l1x, d1;
	if (!p1.getLambda3D(l1y, l1z, l1x, d1)) {
		if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();
		return 0;
	}

	const T t1x = (p2y-p3y);
	const T t1y = (p3x-p2x);
	const T e2 = (l1x*t1x);
	const T e3 = (l1y*t1y);
	const T e = (e2+e3);
	const T pr1 = (p2x*p3y);
	const T pr2 = (p2y*p3x);
	const T pr = (pr1-pr2);
	const T dpr = (d1*pr);
	const T det = (dpr+e);

	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();

	if constexpr (std::is_same<expansion, T>::value) {
		if (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW)) return INT_MAX;
	}
	return sgn(det);
}

inline int orient2dzx_indirect_IEE(const genericPoint& p1, const double& p2x, const double& p2y, const double& p3x, const double& p3y) {
	int ret;
	if ((ret = orient2dzx_indirect_IEE_t<interval_number, interval_number>(p1, p2x, p2y, p3x, p3y)) != 0) return ret;
	return orient2dzx_indirect_IEE_t<bigfloat, bigfloat>(p1, p2x, p2y, p3x, p3y);
}


template<class PT, class T> static inline int orient2dzx_indirect_IIE_t(const genericPoint& p1, const genericPoint& p2, const PT& op3x, const PT& op3y) {
	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundUP();

	std::conditional_t<(std::is_same<expansion, T>::value), expansionPool, char> pool;

	if constexpr (std::is_same<expansion, T>::value) {
		expansion::initPool(&pool);
		feclearexcept(FE_ALL_EXCEPT);
	}

	T l1y, l1z, l1x, d1;
	if (!p1.getLambda3D(l1y, l1z, l1x, d1)) {
		if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();
		return 0;
	}

	T l2y, l2z, l2x, d2;
	if (!p2.getLambda3D(l2y, l2z, l2x, d2)) {
		if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();
		return 0;
	}

	const T a = (d1*l2x);
	const T b = (d2*l1x);
	const T c = (d1*op3y);
	const T e = (d1*l2y);
	const T f = (d2*l1y);
	const T g = (d1*op3x);
	const T ab = (a-b);
	const T cd = (c-l1y);
	const T ef = (e-f);
	const T gh = (g-l1x);
	const T abcd = (ab*cd);
	const T efgh = (ef*gh);
	const T L = (abcd-efgh);

	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();

	if constexpr (std::is_same<expansion, T>::value) {
		if (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW)) return INT_MAX;
	}
	return sgn(L);
}

inline int orient2dzx_indirect_IIE(const genericPoint& p1, const genericPoint& p2, const double& op3x, const double& op3y) {
	int ret;
	if ((ret = orient2dzx_indirect_IIE_t<interval_number, interval_number>(p1, p2, op3x, op3y)) != 0) return ret;
	return orient2dzx_indirect_IIE_t<bigfloat, bigfloat>(p1, p2, op3x, op3y);
}


template<class PT, class T> static inline int orient2dzx_indirect_III_t(const genericPoint& p1, const genericPoint& p2, const genericPoint& p3) {
	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundUP();

	std::conditional_t<(std::is_same<expansion, T>::value), expansionPool, char> pool;

	if constexpr (std::is_same<expansion, T>::value) {
		expansion::initPool(&pool);
		feclearexcept(FE_ALL_EXCEPT);
	}

	T l1y, l1z, l1x, d1;
	if (!p1.getLambda3D(l1y, l1z, l1x, d1)) {
		if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();
		return 0;
	}

	T l2y, l2z, l2x, d2;
	if (!p2.getLambda3D(l2y, l2z, l2x, d2)) {
		if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();
		return 0;
	}

	T l3y, l3z, l3x, d3;
	if (!p3.getLambda3D(l3y, l3z, l3x, d3)) {
		if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();
		return 0;
	}

	const T a = (d1*l2x);
	const T b = (d2*l1x);
	const T c = (d1*l3y);
	const T d = (d3*l1y);
	const T e = (d1*l2y);
	const T f = (d2*l1y);
	const T g = (d1*l3x);
	const T h = (d3*l1x);
	const T ab = (a-b);
	const T cd = (c-d);
	const T ef = (e-f);
	const T gh = (g-h);
	const T abcd = (ab*cd);
	const T efgh = (ef*gh);
	const T L = (abcd-efgh);

	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();

	if constexpr (std::is_same<expansion, T>::value) {
		if (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW)) return INT_MAX;
	}
	return sgn(L);
}

inline int orient2dzx_indirect_III(const genericPoint& p1, const genericPoint& p2, const genericPoint& p3) {
	int ret;
	if ((ret = orient2dzx_indirect_III_t<interval_number, interval_number>(p1, p2, p3)) != 0) return ret;
	return orient2dzx_indirect_III_t<bigfloat, bigfloat>(p1, p2, p3);
}


template<class PT, class T> static inline int orient2d_indirect_IEE_t(const genericPoint& p1, const PT& p2x, const PT& p2y, const PT& p3x, const PT& p3y) {
	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundUP();

	std::conditional_t<(std::is_same<expansion, T>::value), expansionPool, char> pool;

	if constexpr (std::is_same<expansion, T>::value) {
		expansion::initPool(&pool);
		feclearexcept(FE_ALL_EXCEPT);
	}

	T l1x, l1y, d1;
	if (!p1.getLambda2D(l1x, l1y, d1)) {
		if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();
		return 0;
	}

	const T t1x = (p2y-p3y);
	const T t1y = (p3x-p2x);
	const T e2 = (l1x*t1x);
	const T e3 = (l1y*t1y);
	const T e = (e2+e3);
	const T pr1 = (p2x*p3y);
	const T pr2 = (p2y*p3x);
	const T pr = (pr1-pr2);
	const T dpr = (d1*pr);
	const T det = (dpr+e);

	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();

	if constexpr (std::is_same<expansion, T>::value) {
		if (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW)) return INT_MAX;
	}
	return sgn(det);
}

inline int orient2d_indirect_IEE(const genericPoint& p1, const double& p2x, const double& p2y, const double& p3x, const double& p3y) {
	int ret;
	if ((ret = orient2d_indirect_IEE_t<interval_number, interval_number>(p1, p2x, p2y, p3x, p3y)) != 0) return ret;
	return orient2d_indirect_IEE_t<bigfloat, bigfloat>(p1, p2x, p2y, p3x, p3y);
}


template<class PT, class T> static inline int orient2d_indirect_IIE_t(const genericPoint& p1, const genericPoint& p2, const PT& p3x, const PT& p3y) {
	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundUP();

	std::conditional_t<(std::is_same<expansion, T>::value), expansionPool, char> pool;

	if constexpr (std::is_same<expansion, T>::value) {
		expansion::initPool(&pool);
		feclearexcept(FE_ALL_EXCEPT);
	}

	T l1x, l1y, d1;
	if (!p1.getLambda2D(l1x, l1y, d1)) {
		if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();
		return 0;
	}

	T l2x, l2y, d2;
	if (!p2.getLambda2D(l2x, l2y, d2)) {
		if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();
		return 0;
	}

	const T a = (d1*l2x);
	const T b = (d2*l1x);
	const T c = (d1*p3y);
	const T e = (d1*l2y);
	const T f = (d2*l1y);
	const T g = (d1*p3x);
	const T ab = (a-b);
	const T cd = (c-l1y);
	const T ef = (e-f);
	const T gh = (g-l1x);
	const T abcd = (ab*cd);
	const T efgh = (ef*gh);
	const T L = (abcd-efgh);

	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();

	if constexpr (std::is_same<expansion, T>::value) {
		if (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW)) return INT_MAX;
	}
	return sgn(L);
}

inline int orient2d_indirect_IIE(const genericPoint& p1, const genericPoint& p2, const double& p3x, const double& p3y) {
	int ret;
	if ((ret = orient2d_indirect_IIE_t<interval_number, interval_number>(p1, p2, p3x, p3y)) != 0) return ret;
	return orient2d_indirect_IIE_t<bigfloat, bigfloat>(p1, p2, p3x, p3y);
}


template<class PT, class T> static inline int orient2d_indirect_III_t(const genericPoint& p1, const genericPoint& p2, const genericPoint& p3) {
	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundUP();

	std::conditional_t<(std::is_same<expansion, T>::value), expansionPool, char> pool;

	if constexpr (std::is_same<expansion, T>::value) {
		expansion::initPool(&pool);
		feclearexcept(FE_ALL_EXCEPT);
	}

	T l1x, l1y, d1;
	if (!p1.getLambda2D(l1x, l1y, d1)) {
		if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();
		return 0;
	}

	T l2x, l2y, d2;
	if (!p2.getLambda2D(l2x, l2y, d2)) {
		if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();
		return 0;
	}

	T l3x, l3y, d3;
	if (!p3.getLambda2D(l3x, l3y, d3)) {
		if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();
		return 0;
	}

	const T a = (d1*l2x);
	const T b = (d2*l1x);
	const T c = (d1*l3y);
	const T d = (d3*l1y);
	const T e = (d1*l2y);
	const T f = (d2*l1y);
	const T g = (d1*l3x);
	const T h = (d3*l1x);
	const T ab = (a-b);
	const T cd = (c-d);
	const T ef = (e-f);
	const T gh = (g-h);
	const T abcd = (ab*cd);
	const T efgh = (ef*gh);
	const T L = (abcd-efgh);

	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();

	if constexpr (std::is_same<expansion, T>::value) {
		if (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW)) return INT_MAX;
	}
	return sgn(L);
}

inline int orient2d_indirect_III(const genericPoint& p1, const genericPoint& p2, const genericPoint& p3) {
	int ret;
	if ((ret = orient2d_indirect_III_t<interval_number, interval_number>(p1, p2, p3)) != 0) return ret;
	return orient2d_indirect_III_t<bigfloat, bigfloat>(p1, p2, p3);
}


template<class PT, class T> static inline int orient3d_indirect_IEEE_t(const genericPoint& p1, const PT& p2x, const PT& p2y, const PT& p2z, const PT& p3x, const PT& p3y, const PT& p3z, const PT& p4x, const PT& p4y, const PT& p4z) {
	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundUP();

	std::conditional_t<(std::is_same<expansion, T>::value), expansionPool, char> pool;

	if constexpr (std::is_same<expansion, T>::value) {
		expansion::initPool(&pool);
		feclearexcept(FE_ALL_EXCEPT);
	}

	T l1x, l1y, l1z, d1;
	if (!p1.getLambda3D(l1x, l1y, l1z, d1)) {
		if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();
		return 0;
	}

	const T v1x = (p2x-p4x);
	const T v1y = (p2y-p4y);
	const T v1z = (p2z-p4z);
	const T v2x = (p3x-p4x);
	const T v2y = (p3y-p4y);
	const T v2z = (p3z-p4z);
	const T v3x = (l1x-(p4x*d1));
	const T v3y = (l1y-(p4y*d1));
	const T v3z = (l1z-(p4z*d1));
	const T det = (((((v3x*v1z)-(v3z*v1x))*v2y)-(((v3x*v1y)-(v3y*v1x))*v2z))-(((v3y*v1z)-(v3z*v1y))*v2x));

	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();

	if constexpr (std::is_same<expansion, T>::value) {
		if (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW)) return INT_MAX;
	}
	return sgn(det);
}

inline int orient3d_indirect_IEEE(const genericPoint& p1, const double& p2x, const double& p2y, const double& p2z, const double& p3x, const double& p3y, const double& p3z, const double& p4x, const double& p4y, const double& p4z) {
	int ret;
	if ((ret = orient3d_indirect_IEEE_t<interval_number, interval_number>(p1, p2x, p2y, p2z, p3x, p3y, p3z, p4x, p4y, p4z)) != 0) return ret;
	return orient3d_indirect_IEEE_t<bigfloat, bigfloat>(p1, p2x, p2y, p2z, p3x, p3y, p3z, p4x, p4y, p4z);
}


template<class PT, class T> static inline int orient3d_indirect_IIEE_t(const genericPoint& p1, const genericPoint& p2, const PT& p3x, const PT& p3y, const PT& p3z, const PT& p4x, const PT& p4y, const PT& p4z) {
	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundUP();

	std::conditional_t<(std::is_same<expansion, T>::value), expansionPool, char> pool;

	if constexpr (std::is_same<expansion, T>::value) {
		expansion::initPool(&pool);
		feclearexcept(FE_ALL_EXCEPT);
	}

	T l1x, l1y, l1z, d1;
	if (!p1.getLambda3D(l1x, l1y, l1z, d1)) {
		if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();
		return 0;
	}

	T l2x, l2y, l2z, d2;
	if (!p2.getLambda3D(l2x, l2y, l2z, d2)) {
		if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();
		return 0;
	}

	const T v1x = (l2x-(p4x*d2));
	const T v1y = (l2y-(p4y*d2));
	const T v1z = (l2z-(p4z*d2));
	const T v2x = (p3x-p4x);
	const T v2y = (p3y-p4y);
	const T v2z = (p3z-p4z);
	const T v3x = (l1x-(p4x*d1));
	const T v3y = (l1y-(p4y*d1));
	const T v3z = (l1z-(p4z*d1));
	const T det = (((((v3x*v1z)-(v3z*v1x))*v2y)-(((v3x*v1y)-(v3y*v1x))*v2z))-(((v3y*v1z)-(v3z*v1y))*v2x));

	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();

	if constexpr (std::is_same<expansion, T>::value) {
		if (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW)) return INT_MAX;
	}
	return sgn(det);
}

inline int orient3d_indirect_IIEE(const genericPoint& p1, const genericPoint& p2, const double& p3x, const double& p3y, const double& p3z, const double& p4x, const double& p4y, const double& p4z) {
	int ret;
	if ((ret = orient3d_indirect_IIEE_t<interval_number, interval_number>(p1, p2, p3x, p3y, p3z, p4x, p4y, p4z)) != 0) return ret;
	return orient3d_indirect_IIEE_t<bigfloat, bigfloat>(p1, p2, p3x, p3y, p3z, p4x, p4y, p4z);
}


template<class PT, class T> static inline int orient3d_indirect_IIIE_t(const genericPoint& p1, const genericPoint& p2, const genericPoint& p3, const PT& p4x, const PT& p4y, const PT& p4z) {
	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundUP();

	std::conditional_t<(std::is_same<expansion, T>::value), expansionPool, char> pool;

	if constexpr (std::is_same<expansion, T>::value) {
		expansion::initPool(&pool);
		feclearexcept(FE_ALL_EXCEPT);
	}

	T l1x, l1y, l1z, d1;
	if (!p1.getLambda3D(l1x, l1y, l1z, d1)) {
		if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();
		return 0;
	}

	T l2x, l2y, l2z, d2;
	if (!p2.getLambda3D(l2x, l2y, l2z, d2)) {
		if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();
		return 0;
	}

	T l3x, l3y, l3z, d3;
	if (!p3.getLambda3D(l3x, l3y, l3z, d3)) {
		if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();
		return 0;
	}

	const T v1x = (l2x-(p4x*d2));
	const T v1y = (l2y-(p4y*d2));
	const T v1z = (l2z-(p4z*d2));
	const T v2x = (l3x-(p4x*d3));
	const T v2y = (l3y-(p4y*d3));
	const T v2z = (l3z-(p4z*d3));
	const T v3x = (l1x-(p4x*d1));
	const T v3y = (l1y-(p4y*d1));
	const T v3z = (l1z-(p4z*d1));
	const T det = (((((v3x*v1z)-(v3z*v1x))*v2y)-(((v3x*v1y)-(v3y*v1x))*v2z))-(((v3y*v1z)-(v3z*v1y))*v2x));

	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();

	if constexpr (std::is_same<expansion, T>::value) {
		if (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW)) return INT_MAX;
	}
	return sgn(det);
}

inline int orient3d_indirect_IIIE(const genericPoint& p1, const genericPoint& p2, const genericPoint& p3, const double& p4x, const double& p4y, const double& p4z) {
	int ret;
	if ((ret = orient3d_indirect_IIIE_t<interval_number, interval_number>(p1, p2, p3, p4x, p4y, p4z)) != 0) return ret;
	return orient3d_indirect_IIIE_t<bigfloat, bigfloat>(p1, p2, p3, p4x, p4y, p4z);
}


template<class PT, class T> static inline int orient3d_indirect_IIII_t(const genericPoint& p1, const genericPoint& p2, const genericPoint& p3, const genericPoint& p4) {
	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundUP();

	std::conditional_t<(std::is_same<expansion, T>::value), expansionPool, char> pool;

	if constexpr (std::is_same<expansion, T>::value) {
		expansion::initPool(&pool);
		feclearexcept(FE_ALL_EXCEPT);
	}

	T l1x, l1y, l1z, d1;
	if (!p1.getLambda3D(l1x, l1y, l1z, d1)) {
		if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();
		return 0;
	}

	T l2x, l2y, l2z, d2;
	if (!p2.getLambda3D(l2x, l2y, l2z, d2)) {
		if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();
		return 0;
	}

	T l3x, l3y, l3z, d3;
	if (!p3.getLambda3D(l3x, l3y, l3z, d3)) {
		if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();
		return 0;
	}

	T l4x, l4y, l4z, d4;
	if (!p4.getLambda3D(l4x, l4y, l4z, d4)) {
		if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();
		return 0;
	}

	const T v1x = ((l2x*d4)-(l4x*d2));
	const T v1y = ((l2y*d4)-(l4y*d2));
	const T v1z = ((l2z*d4)-(l4z*d2));
	const T v2x = ((l3x*d4)-(l4x*d3));
	const T v2y = ((l3y*d4)-(l4y*d3));
	const T v2z = ((l3z*d4)-(l4z*d3));
	const T v3x = ((l1x*d4)-(l4x*d1));
	const T v3y = ((l1y*d4)-(l4y*d1));
	const T v3z = ((l1z*d4)-(l4z*d1));
	const T det = (((((v3x*v1z)-(v3z*v1x))*v2y)-(((v3x*v1y)-(v3y*v1x))*v2z))-(((v3y*v1z)-(v3z*v1y))*v2x));

	if constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();

	if constexpr (std::is_same<expansion, T>::value) {
		if (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW)) return INT_MAX;
	}
	return sgn(det);
}

inline int orient3d_indirect_IIII(const genericPoint& p1, const genericPoint& p2, const genericPoint& p3, const genericPoint& p4) {
	int ret;
	if ((ret = orient3d_indirect_IIII_t<interval_number, interval_number>(p1, p2, p3, p4)) != 0) return ret;
	return orient3d_indirect_IIII_t<bigfloat, bigfloat>(p1, p2, p3, p4);
}


