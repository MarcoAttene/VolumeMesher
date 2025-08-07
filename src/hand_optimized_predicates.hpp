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

/* Should include incircle and insphere toexpansionObject:: */


#include "implicit_point.h"

#pragma intrinsic(fabs)


inline int orient2d_filtered(double p1x, double p1y, double p2x, double p2y, double p3x, double p3y)
{
	double dl = (p2x - p1x) * (p3y - p1y);
	double dr = (p2y - p1y) * (p3x - p1x);
	double det = dl - dr;
	double eb = 3.3306690738754706e-016 * (fabs(dl) + fabs(dr));
	return ((det >= eb) - (-det >= eb));
}

inline int orient2d_exact(double p1x, double p1y, double p2x, double p2y, double p3x, double p3y)
{
	double acx[2], acy[2], bcx[2], bcy[2], dtl[2], dtr[2], B[4];
	double s[2], t[2], u[4], C1[8], C2[12], D[16];
	int C1l, C2l, Dl;


	acx[1] = (p1x - p3x);
	bcx[1] = (p2x - p3x);
	acy[1] = (p1y - p3y);
	bcy[1] = (p2y - p3y);

	expansionObject::Two_Prod(acx[1], bcy[1], dtl);
	expansionObject::Two_Prod(acy[1], bcx[1], dtr);
	expansionObject::Two_Two_Diff(dtl, dtr, B);

	double dsm = (fabs(dtl[1]) + fabs(dtr[1]));
	double det = expansionObject::To_Double(4, B);
	double eb = 2.2204460492503146e-16 * dsm;
	Dl = ((det >= eb) - (-det >= eb));
	if (Dl) return Dl;

	expansionObject::Two_Diff_Back(p1x, p3x, acx);
	expansionObject::Two_Diff_Back(p2x, p3x, bcx);
	expansionObject::Two_Diff_Back(p1y, p3y, acy);
	expansionObject::Two_Diff_Back(p2y, p3y, bcy);

	if ((acx[0] == 0.0) && (acy[0] == 0.0) && (bcx[0] == 0.0) && (bcy[0] == 0.0)) return ((det > 0) - (det < 0));

	eb = 1.1093356479670487e-31 * dsm + 3.3306690738754706e-16 * fabs(det);
	det += (acx[1] * bcy[0] + bcy[1] * acx[0]) - (acy[1] * bcx[0] + bcx[1] * acy[0]);
	Dl = ((det >= eb) - (-det >= eb));
	if (Dl) return Dl;

	expansionObject::Two_Prod(acx[0], bcy[1], s);
	expansionObject::Two_Prod(acy[0], bcx[1], t);
	expansionObject::Two_Two_Diff(s, t, u);
	C1l = expansionObject::Gen_Sum(4, B, 4, u, C1);

	expansionObject::Two_Prod(acx[1], bcy[0], s);
	expansionObject::Two_Prod(acy[1], bcx[0], t);
	expansionObject::Two_Two_Diff(s, t, u);
	C2l = expansionObject::Gen_Sum(C1l, C1, 4, u, C2);

	expansionObject::Two_Prod(acx[0], bcy[0], s);
	expansionObject::Two_Prod(acy[0], bcx[0], t);
	expansionObject::Two_Two_Diff(s, t, u);
	Dl = expansionObject::Gen_Sum(C2l, C2, 4, u, D);

	det = D[Dl - 1];
	return ((det > 0) - (det < 0));
}

inline int orient2d(double p1x, double p1y, double p2x, double p2y, double p3x, double p3y)
{
	int ret = orient2d_filtered(p1x, p1y, p2x, p2y, p3x, p3y);
	if (ret) return ret;
	return orient2d_exact(p1x, p1y, p2x, p2y, p3x, p3y);
}

inline int orient3d_filtered(double px, double py, double pz, double qx, double qy, double qz, double rx, double ry, double rz, double sx, double sy, double sz)
{
	double fadx, fbdx, fcdx, fady, fbdy, fcdy, fadz, fbdz, fcdz, eb;
	double fbdxcdy, fcdxbdy, fcdxady, fadxcdy, fadxbdy, fbdxady, det;

	fadx = qx - px; fbdx = rx - px; fcdx = sx - px;
	fady = qy - py; fbdy = ry - py; fcdy = sy - py;
	fadz = qz - pz; fbdz = rz - pz; fcdz = sz - pz;

	fbdxcdy = fbdx * fcdy * fadz; fcdxbdy = fcdx * fbdy * fadz;
	fcdxady = fcdx * fady * fbdz; fadxcdy = fadx * fcdy * fbdz;
	fadxbdy = fadx * fbdy * fcdz; fbdxady = fbdx * fady * fcdz;

	det = (fbdxcdy - fcdxbdy) + (fcdxady - fadxcdy) + (fadxbdy - fbdxady);
	eb = 7.7715611723761027e-016 * (fabs(fbdxcdy) + fabs(fcdxbdy) + fabs(fcdxady) + fabs(fadxcdy) + fabs(fadxbdy) + fabs(fbdxady));
	return ((det >= eb) - (-det >= eb));
}

inline int orient3d_bigfloat(const bigfloat px, const bigfloat py, const bigfloat pz, const bigfloat qx, const bigfloat qy, const bigfloat qz, const bigfloat rx, const bigfloat ry, const bigfloat rz, const bigfloat sx, const bigfloat sy, const bigfloat sz)
{
	bigfloat fadx = qx - px; bigfloat fbdx = rx - px; bigfloat fcdx = sx - px;
	bigfloat fady = qy - py; bigfloat fbdy = ry - py; bigfloat fcdy = sy - py;
	bigfloat fadz = qz - pz; bigfloat fbdz = rz - pz; bigfloat fcdz = sz - pz;

	bigfloat fbdxcdy = fbdx * fcdy; bigfloat fcdxbdy = fcdx * fbdy;
	bigfloat fcdxady = fcdx * fady; bigfloat fadxcdy = fadx * fcdy;
	bigfloat fadxbdy = fadx * fbdy; bigfloat fbdxady = fbdx * fady;

	bigfloat det = (fbdxcdy - fcdxbdy) * fadz + (fcdxady - fadxcdy) * fbdz + (fadxbdy - fbdxady) * fcdz;
	return sgn(det);
}

inline int orient3d_exact(double pdx, double pdy, double pdz, double pax, double pay, double paz, double pbx, double pby, double pbz, double pcx, double pcy, double pcz)
{
	double eb, det;
	double adx[2], bdx[2], cdx[2], ady[2], bdy[2], cdy[2], adz[2], bdz[2], cdz[2];
	double bdxcdy[2], cdxbdy[2], cdxady[2], adxcdy[2], adxbdy[2], bdxady[2];
	double bc[4], ca[4], ab[4];
	double adet[8], bdet[8], cdet[8], abdet[16];
	double fin[2][192];
	int wh = 0;
	int alen, blen, clen, finlen;
	int ablen;
	int ri;

	adx[1] = pax - pdx;
	bdx[1] = pbx - pdx;
	cdx[1] = pcx - pdx;
	ady[1] = pay - pdy;
	bdy[1] = pby - pdy;
	cdy[1] = pcy - pdy;
	adz[1] = paz - pdz;
	bdz[1] = pbz - pdz;
	cdz[1] = pcz - pdz;

	expansionObject::Two_Prod(bdx[1], cdy[1], bdxcdy);
	expansionObject::Two_Prod(cdx[1], bdy[1], cdxbdy);
	expansionObject::Two_Two_Diff(bdxcdy, cdxbdy, bc);
	alen = expansionObject::Gen_Scale(4, bc, adz[1], adet);

	expansionObject::Two_Prod(cdx[1], ady[1], cdxady);
	expansionObject::Two_Prod(adx[1], cdy[1], adxcdy);
	expansionObject::Two_Two_Diff(cdxady, adxcdy, ca);
	blen = expansionObject::Gen_Scale(4, ca, bdz[1], bdet);

	expansionObject::Two_Prod(adx[1], bdy[1], adxbdy);
	expansionObject::Two_Prod(bdx[1], ady[1], bdxady);
	expansionObject::Two_Two_Diff(adxbdy, bdxady, ab);
	clen = expansionObject::Gen_Scale(4, ab, cdz[1], cdet);

	ablen = expansionObject::Gen_Sum(alen, adet, blen, bdet, abdet);
	finlen = expansionObject::Gen_Sum(ablen, abdet, clen, cdet, fin[wh]);

	double xx1 = bdxcdy[1] * adz[1]; double xx2 = cdxbdy[1] * adz[1];
	double yy1 = cdxady[1] * bdz[1]; double yy2 = adxcdy[1] * bdz[1];
	double zz1 = adxbdy[1] * cdz[1]; double zz2 = bdxady[1] * cdz[1];
	double pm = fabs(xx1) + fabs(xx2) + fabs(yy1) + fabs(yy2) + fabs(zz1) + fabs(zz2);

	det = expansionObject::To_Double(finlen, fin[wh]);
	eb = 3.3306690738754731e-016 * pm;
	ri = (det >= eb) - (-det >= eb);
	if (ri) return ri;

	expansionObject::Two_Diff_Back(pax, pdx, adx);
	expansionObject::Two_Diff_Back(pbx, pdx, bdx);
	expansionObject::Two_Diff_Back(pcx, pdx, cdx);
	expansionObject::Two_Diff_Back(pay, pdy, ady);
	expansionObject::Two_Diff_Back(pby, pdy, bdy);
	expansionObject::Two_Diff_Back(pcy, pdy, cdy);
	expansionObject::Two_Diff_Back(paz, pdz, adz);
	expansionObject::Two_Diff_Back(pbz, pdz, bdz);
	expansionObject::Two_Diff_Back(pcz, pdz, cdz);

	if ((adx[0] == 0.0) && (bdx[0] == 0.0) && (cdx[0] == 0.0) &&
		(ady[0] == 0.0) && (bdy[0] == 0.0) && (cdy[0] == 0.0) &&
		(adz[0] == 0.0) && (bdz[0] == 0.0) && (cdz[0] == 0.0)) return (det > 0) - (det < 0);

	eb = 3.2047474274603644e-031 * pm + 1.1102230246251565e-016 * fabs(det);
	det += (adz[1] * ((bdx[1] * cdy[0] + cdy[1] * bdx[0])
		- (bdy[1] * cdx[0] + cdx[1] * bdy[0]))
		+ adz[0] * (bdx[1] * cdy[1] - bdy[1] * cdx[1]))
		+ (bdz[1] * ((cdx[1] * ady[0] + ady[1] * cdx[0])
			- (cdy[1] * adx[0] + adx[1] * cdy[0]))
			+ bdz[0] * (cdx[1] * ady[1] - cdy[1] * adx[1]))
		+ (cdz[1] * ((adx[1] * bdy[0] + bdy[1] * adx[0])
			- (ady[1] * bdx[0] + bdx[1] * ady[0]))
			+ cdz[0] * (adx[1] * bdy[1] - ady[1] * bdx[1]));
	ri = (det >= eb) - (-det >= eb);
	if (ri) return ri;

	// Filters did not work. Compute exactly...
	return orient3d_bigfloat(pdx, pdy, pdz, pax, pay, paz, pbx, pby, pbz, pcx, pcy, pcz);
}

inline int orient3d(double px, double py, double pz, double qx, double qy, double qz, double rx, double ry, double rz, double sx, double sy, double sz)
{
	int ret;
	ret = orient3d_filtered(px, py, pz, qx, qy, qz, rx, ry, rz, sx, sy, sz);
	if (ret) return ret;
	return orient3d_exact(px, py, pz, qx, qy, qz, rx, ry, rz, sx, sy, sz);
}
