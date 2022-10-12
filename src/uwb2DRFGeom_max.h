#pragma once

using namespace gismo;

template<class TT> void defineBCs_NS(gsBoundaryConditions<TT>& bcInfo, TT velocity_x, TT velocity_y, TT velocity_blade)
{

    std::ostringstream strs_vel_x;    strs_vel_x << velocity_x;        std::string str_x = strs_vel_x.str();
    std::ostringstream strs_velblade; strs_velblade << velocity_blade; std::string str2 = strs_velblade.str();
    std::ostringstream strs_vel_y;    strs_vel_y << velocity_y;        std::string str_y = strs_vel_y.str();
    
    gsFunctionExpr<TT> Uin(str_x, str_y, 2);
    gsFunctionExpr<TT> Ublade("0",str2, 2);
    //gsFunctionExpr<TT> Uwall("0", "0", 2);	

        //bcInfo.addCondition(0, boundary::south, condition_type::dirichlet, Uwall, 0);
        //bcInfo.addCondition(0, boundary::north, condition_type::dirichlet, Uwall, 0);
        bcInfo.addCondition(1, boundary::south, condition_type::dirichlet, Ublade, 0);
        bcInfo.addCondition(1, boundary::north, condition_type::dirichlet, Ublade, 0);
        //bcInfo.addCondition(2, boundary::south, condition_type::dirichlet, Ublade, 0);
        //bcInfo.addCondition(2, boundary::north, condition_type::dirichlet, Ublade, 0);
        bcInfo.addCondition(0, boundary::west, condition_type::dirichlet, Uin, 0);
        //---------2D periodic conditions------
        //bcInfo.setIdentityMatrix(2);
        //bcInfo.addPeriodic(0, boundary::north, 0, boundary::south, 2);
        //bcInfo.addPeriodic(2, boundary::north, 2, boundary::south, 2);
}

template<class TT> void defineBCs_TM(gsBoundaryConditions<TT>& bcInfoTurb, TT kIn, TT kWall, TT oIn, TT oBlade)
{
    // Boundary conditions
    std::ostringstream str_kIn;    str_kIn << kIn;       std::string s_kIn = str_kIn.str();       gsFunctionExpr<TT> Kin(s_kIn, 2);
    std::ostringstream str_oIn;    str_oIn << oIn;       std::string s_oIn = str_oIn.str();       gsFunctionExpr<TT> Oin(s_oIn, 2);
    std::ostringstream str_kWall;  str_kWall << kWall;   std::string s_kWall = str_kWall.str();   gsFunctionExpr<TT> Kwall(s_kWall, 2);
//    std::ostringstream str_oWall;  str_oWall << oWall;   std::string s_oWall = str_oWall.str();   gsFunctionExpr<TT> Owall(s_oWall, 2);
    std::ostringstream str_oBlade; str_oBlade << oBlade; std::string s_oBlade = str_oBlade.str(); gsFunctionExpr<TT> Oblade(s_oBlade, 2);
        bcInfoTurb.addCondition(1, boundary::south, condition_type::dirichlet, Kwall, 0);
        bcInfoTurb.addCondition(1, boundary::north, condition_type::dirichlet, Kwall, 0);
        bcInfoTurb.addCondition(0, boundary::west, condition_type::dirichlet, Kin, 0);

        bcInfoTurb.addCondition(1, boundary::south, condition_type::dirichlet, Oblade, 1);
        bcInfoTurb.addCondition(1, boundary::north, condition_type::dirichlet, Oblade, 1);
        bcInfoTurb.addCondition(0, boundary::west, condition_type::dirichlet, Oin, 1);
}

template<class TT> void refineBasis(gsMultiBasis<TT>& tbasis, int numRefineUniform, int numRefineLocal)
{
        gsMatrix<TT> box_v(2,2);
        gsMatrix<TT> box_u(2,2);

        //uniform refine
        for (int i = 0; i < numRefineUniform; i++)
        {
            tbasis.uniformRefine();
        }

        const gsTensorBSplineBasis<2, TT>*  basis = dynamic_cast<const gsTensorBSplineBasis<2, TT>*>(&(tbasis.basis(0))); //basis of the first patch
        TT uRefineKnot = basis->knot(1,basis->degree(1)+1); // first non-zero knot in direction v
        //gsInfo << "uRefineKnot: " << uRefineKnot << "\n";

        //refine in v, near uper wall
        for (int i = 0; i < numRefineLocal; i++)
        {
            box_v << 0.0, 0.0, 1.0 - (uRefineKnot / math::pow(2, i)), 1.0;

            tbasis.refine(0, box_v);
            tbasis.refine(1, box_v);
            tbasis.refine(2, box_v);
        }


        //refine in v, near bottom wall
        for (int i = 0; i < numRefineLocal; i++)
        {
            box_v << 0.0, 0.0, 0.0, uRefineKnot / math::pow(2, i);

            tbasis.refine(0, box_v);
            tbasis.refine(1, box_v);
            tbasis.refine(2, box_v);
        }

        //refine in u, before blade
        for (int i = 0; i < numRefineLocal; i++)
        {
            box_u << 1.0 - uRefineKnot / math::pow(2, i), 1.0, 0.0, 0.0;

            tbasis.refine(0, box_u);
        }

        //refine in u, blade
        for (int i = 0; i < numRefineLocal; i++)
        {
            box_u << 0.0, uRefineKnot / math::pow(2, i), 0.0 , 0.0;

            tbasis.refine(1, box_u);
        }
}

template<class TT> void refineBasisUniformZones(gsMultiBasis<TT>& tbasis, int numRefineUniformLocal_v,
                                               int numUniformKnot_v)
{
    gsMatrix<TT> box_v_bottom(2,2);
    gsMatrix<TT> box_v_up(2,2);

    for (int i = 0; i < numRefineUniformLocal_v; i++)
    {
         // in each refinement step take the first numKnot_v non_zero knots
         const gsTensorBSplineBasis<2, TT>*  basis_0 = dynamic_cast<const gsTensorBSplineBasis<2, TT>*>(&(tbasis.basis(0))); //basis of the first patch
         gsVector<TT> uRefineKnot_v_start((i+1)*(numUniformKnot_v)+1);
         gsVector<TT> uRefineKnot_v_end((i+1)*(numUniformKnot_v)+1);
         uRefineKnot_v_start.setZero((i+1)*(numUniformKnot_v)+1);
         uRefineKnot_v_end.setZero((i+1)*(numUniformKnot_v)+1);
         int sizeKnots_v = basis_0->knots(1).size()-1;
         for (int k = 0; k < (i+1)*(numUniformKnot_v)+1; k++) //gsVector of first numKnots knots
         {
             uRefineKnot_v_start(k) = basis_0->knot(1,basis_0->degree(1)+k); // first numKnot knots and 0 knot in direction v - all patches have same knots in v direction
              uRefineKnot_v_end(k) = basis_0->knot(1,sizeKnots_v - (basis_0->degree(1)+k)); //last numKnot knots and 1 knot in direction v
         }

         for (int j = 0; j < (i+1)*(numUniformKnot_v); j++)
         {
            box_v_bottom << 0.0, 0.0, uRefineKnot_v_start(j), uRefineKnot_v_start(j + 1);
            box_v_up << 0.0, 0.0, uRefineKnot_v_end(j + 1), uRefineKnot_v_end(j);

            tbasis.refine(0, box_v_up);
            tbasis.refine(1, box_v_up);
            tbasis.refine(2, box_v_up);
            tbasis.refine(0, box_v_bottom);
            tbasis.refine(1, box_v_bottom);
            tbasis.refine(2, box_v_bottom);
         }
    }
}

template<class TT> void refineBasisZones(gsMultiBasis<TT>& tbasis,
                                        int numRefineLocal_v, int numKnot_v,
                                        int numRefineLocal_u0, int numRefineLocal_u1_s,int numRefineLocal_u1_e,int numRefineLocal_u2, int numKnot_u0, int numKnot_u1_s,  int numKnot_u1_e, int numKnot_u2)
{
        gsMatrix<TT> box_v_bottom(2,2);
        gsMatrix<TT> box_v_up(2,2);
        gsMatrix<TT> box_u0(2,2);
        gsMatrix<TT> box_u1_s(2,2);
        gsMatrix<TT> box_u1_e(2,2);
        gsMatrix<TT> box_u2(2,2);

        //knots in v direction assumption !!!SAME FOR ALL PATCHES
        gsVector<TT> uRefineKnot_v_start(numKnot_v+1);
        gsVector<TT> uRefineKnot_v_end(numKnot_v+1);
        gsVector<TT> uRefineKnot_u0(numKnot_u0+1);
        gsVector<TT> uRefineKnot_u1_s(numKnot_u1_s+1);
        gsVector<TT> uRefineKnot_u1_e(numKnot_u1_e+1);
        gsVector<TT> uRefineKnot_u2(numKnot_u2+1);

        //refine in v, near bottom and upper wall
        for (int i = 0; i < numRefineLocal_v; i++)
        {
             // in each refinement step take the first numKnot_v non_zero knots
             const gsTensorBSplineBasis<2, TT>*  basis_0 = dynamic_cast<const gsTensorBSplineBasis<2, TT>*>(&(tbasis.basis(0))); //basis of the first patch
             uRefineKnot_v_start.setZero(numKnot_v+1);
             uRefineKnot_v_end.setZero(numKnot_v+1);
             int sizeKnots_v = basis_0->knots(1).size()-1;
             for (int k = 0; k < numKnot_v+1; k++) //gsVector of first numKnots knots
             {
                 uRefineKnot_v_start(k) = basis_0->knot(1,basis_0->degree(1)+k); // first numKnot knots and 0 knot in direction v - all patches have same knots in v direction
                  uRefineKnot_v_end(k) = basis_0->knot(1,sizeKnots_v - (basis_0->degree(1)+k)); //last numKnot knots and 1 knot in direction v
             }

             for (int j = 0; j < numKnot_v; j++)
             {
                box_v_bottom << 0.0, 0.0, uRefineKnot_v_start(j), uRefineKnot_v_start(j + 1);
                box_v_up << 0.0, 0.0, uRefineKnot_v_end(j + 1), uRefineKnot_v_end(j);

                tbasis.refine(0, box_v_up);
                tbasis.refine(1, box_v_up);
                tbasis.refine(2, box_v_up);
                tbasis.refine(0, box_v_bottom);
                tbasis.refine(1, box_v_bottom);
                tbasis.refine(2, box_v_bottom);
             }
        }

        //first patch, refine in u -> before blade
        for (int i = 0; i < numRefineLocal_u0; i++)
        {
             // in each refinement step take the first numKnot_v non_zero knots
             const gsTensorBSplineBasis<2, TT>*  basis_u0 = dynamic_cast<const gsTensorBSplineBasis<2, TT>*>(&(tbasis.basis(0))); //basis of the first patch
             uRefineKnot_u0.setZero(numKnot_u0+1);
             int sizeKnots = basis_u0->knots(0).size()-1;

             for (int k = 0; k < numKnot_u0+1; k++) //gsVector of first numKnots knots
             {
                 uRefineKnot_u0(k) = basis_u0->knot(0,sizeKnots -(basis_u0->degree(0)+k)); // first numKnot knots and 0 knot in direction u first patch      }
             }

             for (int j = 0; j < numKnot_u0; j++)
             {
                 box_u0 <<  uRefineKnot_u0(j + 1), uRefineKnot_u0(j), 0.0, 0.0;
                 tbasis.refine(0, box_u0);
             }
        }

        //second patch, refine in u -> leading edge of blade
        for (int i = 0; i < numRefineLocal_u1_s; i++)
        {
             // in each refinement step take the first numKnot_v non_zero knots
             const gsTensorBSplineBasis<2, TT>*  basis_u1_s = dynamic_cast<const gsTensorBSplineBasis<2, TT>*>(&(tbasis.basis(1))); //basis of the first patch

             uRefineKnot_u1_s.setZero(numKnot_u1_s+1);
             for (int k = 0; k < numKnot_u1_s+1; k++) //gsVector of first numKnots knots
             {
                 uRefineKnot_u1_s(k) = basis_u1_s->knot(0,basis_u1_s->degree(0)+k); // first numKnot knots and 0 knot in direction v - all patches have same knots in v direction        }
             }

             for (int j = 0; j < numKnot_u1_s; j++)
             {
                box_u1_s << uRefineKnot_u1_s(j), uRefineKnot_u1_s(j + 1), 0.0, 0.0;
                tbasis.refine(1, box_u1_s);
             }
        }

        //second patch, refine in u -> end
        for (int i = 0; i < numRefineLocal_u1_e; i++)
        {
             // in each refinement step take the first numKnot_v non_zero knots
             const gsTensorBSplineBasis<2, TT>*  basis_u1_e = dynamic_cast<const gsTensorBSplineBasis<2, TT>*>(&(tbasis.basis(1))); //basis of the first patch
             uRefineKnot_u1_e.setZero(numKnot_u1_e+1);
             int sizeKnots_1 = basis_u1_e->knots(0).size()-1;

             for (int k = 0; k < numKnot_u1_e+1; k++) //gsVector of first numKnots knots
             {
                 uRefineKnot_u1_e(k) = basis_u1_e->knot(0,sizeKnots_1 -(basis_u1_e->degree(0)+k)); // first numKnot knots and 0 knot in direction u first patch      }
             }

             for (int j = 0; j < numKnot_u1_e; j++)
             {
                 box_u1_e <<  uRefineKnot_u1_e(j + 1), uRefineKnot_u1_e(j), 0.0, 0.0;
                 tbasis.refine(1, box_u1_e);
             }
        }

        //last patch, refine in u -> start
        for (int i = 0; i < numRefineLocal_u2; i++)
        {
             // in each refinement step take the first numKnot_v non_zero knots
             const gsTensorBSplineBasis<2, TT>*  basis_u2 = dynamic_cast<const gsTensorBSplineBasis<2, TT>*>(&(tbasis.basis(2))); //basis of the first patch

             uRefineKnot_u2.setZero(numKnot_u2+1);
             for (int k = 0; k < numKnot_u2+1; k++) //gsVector of first numKnots knots
             {
                 uRefineKnot_u2(k) = basis_u2->knot(0,basis_u2->degree(0)+k); // first numKnot knots and 0 knot in direction u - all patches have same knots in v direction        }
             }

             for (int j = 0; j < numKnot_u2; j++)
             {
                box_u2 << uRefineKnot_u2(j), uRefineKnot_u2(j + 1), 0.0, 0.0;
                tbasis.refine(2, box_u2);
             }
        }

}

template<class TT>
gsMultiPatch<TT> BSplineProfile2DBetweenPatch(int const & index,
                                       TT const & length_x1,
                                       TT const & length_x2,
                                       TT const & pitch,
                                       TT const & camberX,
                                       TT const & camberY,
                                       TT const & leadingAngle,
                                       TT const & trailingAngle,
                                       TT const & thicknessX,
                                       TT const & thicknessY,
                                       TT const & endingOffset,
                                       TT const & outputAngle,
                                       TT const & radius,
                                       TT const & chordLength,
                                       TT const & Angle,
                                       TT const & rotationCenterX,
                                       TT const & rotationCenterY,
                                       TT const & uniformity_param)
{

    //----------------set parameters for blade profile----------------
    bool plot = false;
    int num_samples = 30;
    gsVector<TT> vec(2);
    //gsInfo << pitch << "\n ";
    vec (0) = rotationCenterX;
    vec (1) = rotationCenterY;
    gsMatrix<TT> mat(2,2);
    mat << chordLength * math::cos(Angle), chordLength * math::sin(Angle),
           chordLength * math::sin(Angle), -chordLength * math::cos(Angle);
    //gsInfo << vec << "\n \n";
    //gsInfo << mat << "\n \n";
    gsBSpline<TT> suction_side_curve;
    gsBSpline<TT> pressure_side_curve;
    gsBSpline<TT> suction_side_curve_transf;
    gsBSpline<TT> pressure_side_curve_transf;
    gsKnotVector<TT> kvfit(0, 1, 4, 4);
    unsigned num_cpblade = 8;
    BladeProfile<TT> * pBladeProfile = 0;
    //unsigned degree = 3;
    //---------------compute blade profile for given parameters----------------------
    pBladeProfile = new BladeProfile<TT>(camberX, camberY, leadingAngle, trailingAngle, thicknessX,
                                       thicknessY, endingOffset, outputAngle, radius, chordLength,
                                       Angle, rotationCenterX, rotationCenterY, 0.0);
    pBladeProfile->compute2D(suction_side_curve, pressure_side_curve, kvfit, num_samples, 3, 2);
     //---------------transform given profile----------------------
    //gsInfo << suction_side_curve.coefs();
    //gsInfo << pressure_side_curve.coefs();

    suction_side_curve.translate(-vec);
    pressure_side_curve.translate(-vec);
    pBladeProfile->setSuctionSide(suction_side_curve);
    pBladeProfile->setPressureSide(pressure_side_curve);
    pressure_side_curve_transf = pBladeProfile->getPressureSide();
    suction_side_curve_transf = pBladeProfile->getSuctionSide();
    pressure_side_curve_transf.linearTransform(mat);
    suction_side_curve_transf.linearTransform(mat);

    //gsInfo << suction_side_curve_transf.coefs();
    //gsInfo << pressure_side_curve_transf.coefs();

    pBladeProfile->setPressureSide(pressure_side_curve_transf);
    pBladeProfile->setSuctionSide(suction_side_curve_transf);
    gsBSpline < TT > bs = pBladeProfile -> getPressureSide ();
    gsBSpline < TT > bp = pBladeProfile -> getSuctionSide ();
    gsMatrix < TT > cp_bp (num_cpblade, 2);
    gsMatrix < TT > cp_bp_pom (num_cpblade, 2);
    gsMatrix < TT > cp_bs (num_cpblade, 2);

    //---------------set parameters for boundary of patches-----------------------
    //real_t width_y2 = /*bp.coef(bp.coefsSize()-1, 1)/2.0*/ + pitch/2.0;

    //control points of pressure side
    for(unsigned i = 0; i < num_cpblade; i++){
        cp_bp(i, 0) = bp.coef(i, 0) ;
        cp_bp(i, 1) = bp.coef(i, 1) + pitch;
        cp_bp_pom(i, 0) = bp.coef(i, 0) ;
        cp_bp_pom(i, 1) = bp.coef(i, 1) ;
    }
    //control points of suction side
    for(unsigned i = 0; i < num_cpblade; i++){
        cp_bs(i, 0) = bs.coef(i, 0);
        cp_bs(i, 1) = bs.coef(i, 1);
    }

    /*gsInfo << "\n pressure side \n";
    gsInfo << cp_bp;
    gsInfo << "\n pressure side pom \n";
    gsInfo << cp_bp_pom;
    gsInfo << "\n suction side \n";
    gsInfo << cp_bs;
    */

    gsMultiPatch<TT> mp;
    // compute discrete coons patch to optimize
    gsMatrix<TT> coef_patchAll (56, 2);
    coef_patchAll.setZero(56, 2);
    gsMatrix<TT> a_cp(14,2);
    gsMatrix<TT> b_cp(14,2);
    gsMatrix<TT> c_cp(4,2);
    gsMatrix<TT> d_cp(4,2);

    TT ystart_coor = -((cp_bs(0,1)*cp_bs(7,0) - cp_bs(0,0)*cp_bs(7,1) - cp_bs(0,1)*length_x1+ cp_bs(7,1)*length_x1)/(
                          cp_bs(0,0) - cp_bs(7,0)));

    gsMatrix<TT> coef_patchStart (4, 2);
    coef_patchStart << length_x1, ystart_coor,
                         3.0*length_x1/4.0 + cp_bs(0,0)/4.0, 3.0*ystart_coor/4.0 + cp_bs(0,1)/4.0,
                         length_x1/4.0 + 3.0*cp_bs(0,0)/4.0, ystart_coor/4.0 + 3.0*cp_bs(0,1)/4.0,
                         cp_bs(0,0), cp_bs(0,1);

     for (unsigned int i = 0; i<4;i++)
     {
        a_cp(i,0)=coef_patchStart(i,0);
        a_cp(i,1)=coef_patchStart(i,1);
        b_cp(i,0)=coef_patchStart(i,0);
        b_cp(i,1)=coef_patchStart(i,1)+pitch;
     }

    for (unsigned int i = 4; i<10;i++)
    {
        a_cp(i,0)=cp_bs(i-3,0);
        a_cp(i,1)=cp_bs(i-3,1);
        b_cp(i,0)=cp_bp(i-3,0);
        b_cp(i,1)=cp_bp(i-3,1);
    }

    TT yend_coor = -((cp_bs(0,1)*cp_bs(7,0) - cp_bs(0,0)*cp_bs(7,1) - cp_bs(0,1)*length_x2 + cp_bs(7,1)*length_x2)/(
                          cp_bs(0,0) - cp_bs(7,0)));
    gsMatrix<TT> coef_patchEnd (4, 2);
    coef_patchEnd << cp_bs(7,0),cp_bs(7,1),
                     length_x2/4.0 + 3.0*cp_bs(7,0)/4.0, yend_coor/4.0 + 3.0*cp_bs(7,1)/4.0,
                     3.0*length_x2/4.0 + cp_bs(7,0)/4.0, 3.0*yend_coor/4.0 + cp_bs(7,1)/4.0,
                     length_x2, yend_coor;
    for (unsigned int i = 10; i<14;i++)
    {
        a_cp(i,0)=coef_patchEnd(i-10,0);
        a_cp(i,1)=coef_patchEnd(i-10,1);
        b_cp(i,0)=coef_patchEnd(i-10,0);
        b_cp(i,1)=coef_patchEnd(i-10,1)+pitch;
    }
        c_cp << length_x1,a_cp(0,1),
                length_x1,1.0*b_cp(0,1)/4.0 + 3.0*a_cp(0,1)/4.0,
                length_x1,3.0*b_cp(0,1)/4.0 + 1.0*a_cp(0,1)/4.0,
                length_x1, b_cp(0,1);
        d_cp << length_x2,a_cp(13,1),
                length_x2,1.0*b_cp(13,1)/4.0 + 3.0*a_cp(13,1)/4.0,
                length_x2,3.0*b_cp(13,1)/4.0 + 1.0*a_cp(13,1)/4.0,
                length_x2, b_cp(13,1);
    gsKnotVector<TT> kv_uu(0, 1, 4, 4);
    gsKnotVector<TT> kv_vv(0, 1, 0, 4);
    kv_uu.insert(0.2/3.0,3);
    kv_uu.insert(1-(0.2/3.0),3);
    gsInfo<< kv_uu;
    gsMultiPatch<TT> * boundaries4 = new gsMultiPatch<TT>;
    boundaries4->addPatch(gsBSpline<TT>( kv_vv, c_cp));
    boundaries4->addPatch(gsBSpline<TT>( kv_uu, a_cp));
    boundaries4->addPatch(gsBSpline<TT>( kv_vv, d_cp));
    boundaries4->addPatch(gsBSpline<TT>( kv_uu, b_cp));
    gsCoonsPatch<TT> patchAll = coonsPatch(*boundaries4);
    patchAll.compute();
    mp.addPatch(patchAll.result());

   //=================================optimization===========================================

    gsVector<TT> area_vec(7);
    area_vec.setZero(7);
    area_vec << 0.1, 0.25, 0.5, 0.5, 0.5, 0.75, 1.0;//1,1,1,1,0.75,0.75,0.75;

    TT orthogonality = 0.0;
    TT skewness = 0.0;
    TT eccentricity = 0.0;
    TT intersection = 0.0;
    TT uniformity = uniformity_param;
    TT area = area_vec(index);
    TT length = 0.0;
    TT epsilon = 1e-7;

    gsQualityMeasure<TT> optimization(mp.patch(0));
    optimization.optimize(orthogonality, skewness,
                          eccentricity, uniformity,
                          length, area,
                          intersection, epsilon);

    if(plot)
    {
        gsFileData<TT> fileData;
        fileData << mp.patch(0);
        std::string out;
        out = "optimize_blade"+ util::to_string(index) +".xml";
        fileData.dump(out);
        gsMesh<TT> mesh;
        makeMesh(mp.patch(0).basis(), mesh, 10);
        mp.patch(0).evaluateMesh(mesh);
        out = "optimize_bladeMesh" +  util::to_string(index) ;
        gsWriteParaview(mesh, out);
        gsMesh<TT> mesh2;
        mp.patch(0).controlNet(mesh2);
        out = "optimize_bladeControlNet" +  util::to_string(index);
        gsWriteParaview(mesh2,out);
    }

    //=================================divide gsGeometry into three patches===========================================

    //initial data for patches
    gsMatrix<TT> coefs = mp.patch(0).coefs();
    gsMultiPatch<TT> mpFinal;
    gsKnotVector<TT> kv_u(0, 1, 0, 4);
    gsKnotVector<TT> kv_v(0, 1, 0, 4);
    gsTensorBSplineBasis<2, TT> basis(kv_u, kv_v);
    gsTensorBSplineBasis<2, TT> basis_blade(kvfit, kv_v);

    //--------------------------------patch 0-------------------------------------------
    gsMatrix<TT> coef_patch0(16, 2);
    coef_patch0.setZero(16,2);

    for(int i = 0; i < 4; i++)
    {
       for(int j = 0; j < 4; j++)
       {
          coef_patch0(i*4+j,0) = coefs(i*4+j+i*10,0);
          coef_patch0(i*4+j,1) = coefs(i*4+j+i*10,1);
       }
    }
    gsTensorBSpline<2, TT> patch0(basis, coef_patch0);
    for(TT knot = 0.2; knot < 1.0; knot += 0.2)
    {
        patch0.insertKnot(knot,0);
        patch0.insertKnot(knot,1);
    }

    //--------------------------------patch 1-------------------------------------------
    gsMatrix<TT> coef_patch1 (num_cpblade*4, 2);
    coef_patch1.setZero(num_cpblade*4,2);

    for(int i = 0; i < 4; i++)
    {
       for(unsigned j = 0; j < num_cpblade; j++)
       {
          coef_patch1(i*num_cpblade+j,0) = coefs(i*4+j+3+i*10,0);
          coef_patch1(i*num_cpblade+j,1) = coefs(i*4+j+3+i*10,1);
       }
    }
    //gsInfo << "\n pressure side \n";
    //gsInfo << coef_patch1<< "\n";
    gsTensorBSpline<2, TT> patch1(basis_blade, coef_patch1);
    for(TT knot = 0.2; knot < 1.0; knot += 0.2)
    {
        patch1.insertKnot(knot,1);
    }

    //--------------------------------patch 2-------------------------------------------
    gsMatrix<TT> coef_patch2(16, 2);
    coef_patch2.setZero(16,2);

    for(int i = 0; i < 4; i++)
    {
       for(int j = 0; j < 4; j++)
       {
          coef_patch2(i*4+j,0) = coefs(i*4+j+10+i*10,0);
          coef_patch2(i*4+j,1) = coefs(i*4+j+10+i*10,1);
       }
    }
    gsTensorBSpline<2, TT> patch2(basis, coef_patch2);
    for(TT knot = 0.2; knot < 1.0; knot += 0.2)
    {
        patch2.insertKnot(knot,0);
        patch2.insertKnot(knot,1);
    }

    mpFinal.addPatch(patch0);
    mpFinal.addPatch(patch1);
    mpFinal.addPatch(patch2);

    mpFinal.addInterface(0, boundary::east, 1, boundary::west);
    mpFinal.addInterface(1, boundary::east, 2, boundary::west);
    mpFinal.addAutoBoundaries();

    return mpFinal;
}

template<class TT>
void solvePoissonEquation(gsMultiPatch<TT> patches, uwbTMSolverKOmega<TT>& turbSolver_unsteady, std::string profile, bool geometry, int plot_pts)
{
    bool plotMeshes = false;

    int numRefinePoisson;
    if (geometry)
        numRefinePoisson = 4;
    else
        numRefinePoisson = 3;

    int numRefineUniformLocal_vPoisson = 0;
    int numUniformKnot_vPoisson = 1;

    int numRefineLocal_vPoisson = 0;//4;
    int numKnot_vPoisson = 2;
    int numRefineLocal_u0Poisson = 0;//4;
    int numKnot_u0Poisson = 1;
    int numRefineLocal_u1_sPoisson = 0;//4;
    int numKnot_u1_sPoisson = 1;
    int numRefineLocal_u1_ePoisson = 0;
    int numKnot_u1_ePoisson = 1;
    int numRefineLocal_u2Poisson = 0;
    int numKnot_u2Poisson = 1;

    int numRefineLocalFirstKnot_vPoisson = 0;//5; //after refinement one more time from the first knot (refined net) in v

    gsMultiBasis<> tbasisPoisson(patches); // basis for RANS equations
    for (int i = 0; i < numRefinePoisson; ++i)
        tbasisPoisson.uniformRefine();
    refineBasisUniformZones(tbasisPoisson, numRefineUniformLocal_vPoisson, numUniformKnot_vPoisson);
    gsMatrix<> box_u0Poisson(2, 2);
    box_u0Poisson << 0, 1, 0, 0;
    tbasisPoisson.refine(1, box_u0Poisson);
    refineBasisZones(tbasisPoisson, numRefineLocal_vPoisson, math::pow(2, numRefineUniformLocal_vPoisson) * numKnot_vPoisson, numRefineLocal_u0Poisson, numRefineLocal_u1_sPoisson, numRefineLocal_u1_ePoisson, numRefineLocal_u2Poisson, numKnot_u0Poisson, numKnot_u1_sPoisson, numKnot_u1_ePoisson, numKnot_u2Poisson); //for geometry_between = true
                                                                                                                                                                                                                                       //one more local refinement in v from 1st knot (in already refined mesh)
    refineBasisZones(tbasisPoisson, numRefineLocalFirstKnot_vPoisson, 1, 0, 0, 0, 0, 0, 0, 0, 0);
    if (plotMeshes)
    {
        gsMesh<> mesh;
        makeMesh(tbasisPoisson.at(0), mesh, 10);
        patches.patch(0).evaluateMesh(mesh);
        gsWriteParaview(mesh, "meshPoissonPatch0" + profile);
        gsMesh<> mesh1;
        makeMesh(tbasisPoisson.at(1), mesh1, 10);
        patches.patch(1).evaluateMesh(mesh1);
        gsWriteParaview(mesh1, "meshPoissonPatch1" + profile);
        gsMesh<> mesh2;
        makeMesh(tbasisPoisson.at(2), mesh2, 10);
        patches.patch(2).evaluateMesh(mesh2);
        gsWriteParaview(mesh2, "meshPoissonPatch2" + profile);
    }

    gsFunctionExpr<real_t> fw("1", 2);
    gsFunctionExpr<real_t> gw("0", 2);
    gsFunctionExpr<real_t> wallw("0.0", 2);
    gsBoundaryConditions<real_t> bcInfow;
    bcInfow.addCondition(0, boundary::north, condition_type::neumann, gw);
    bcInfow.addCondition(0, boundary::south, condition_type::neumann, gw);
    bcInfow.addCondition(0, boundary::west, condition_type::neumann, gw);
    bcInfow.addCondition(2, boundary::north, condition_type::neumann, gw);
    bcInfow.addCondition(2, boundary::south, condition_type::neumann, gw);
    bcInfow.addCondition(2, boundary::east, condition_type::neumann, gw);
    bcInfow.addCondition(1, boundary::north, condition_type::dirichlet, wallw);
    bcInfow.addCondition(1, boundary::south, condition_type::dirichlet, wallw);

    gsInfo << "\nSolving Poisson equation.\n";

    turbSolver_unsteady.setPoissonSolution(patches, tbasisPoisson, bcInfow, fw, true, plot_pts);

    gsInfo << "Poisson equation resolved.\n\n";
    //turbSolver_unsteady.plotWallDistance("profile_wall_distance" + profile, plot_pts);
}



/*template<class TT>
TT computeWallDistance(gsMultiBasis<TT>& tbasis, gsMultiPatch<TT>&
patches, gsVector<int> distancePatches, std::vector<boxSide>
distanceSides, TT viscosity, TT reynoldsNumber, TT uFreeStream, unsigned npts)
{
    gsVector<TT> minYPlusOverSides;
    minYPlusOverSides.setZero(distancePatches.rows());

    TT skinFrCoeff = math::pow(2.0*math::log10(reynoldsNumber) - 0.65,
-2.3);
    TT tau_w = 0.5 * skinFrCoeff * math::pow(uFreeStream, 2);
    TT u_tau = math::sqrt(tau_w);

    for (int numSide = 0; numSide < distancePatches.rows(); numSide++)
    {
        int patchIndex = distancePatches[numSide];
        boxSide side = distanceSides[numSide];

        const gsTensorBSplineBasis<2, TT>*  basis = dynamic_cast<const
gsTensorBSplineBasis<2, TT>*>(&tbasis.piece(patchIndex));
        typename gsGeometry<TT>::uPtr geomBoundary =
patches.patch(patchIndex).boundary(side);

        gsMatrix<TT> ab = geomBoundary->support();
        gsVector<TT> a = ab.col(0);
        gsVector<TT> b = ab.col(1);

        gsVector<unsigned> np = uniformSampleCount(a, b, npts);
        gsMatrix<TT> paramBoundaryNodes = gsPointGrid(a, b, np);
        int numOfPoints = paramBoundaryNodes.cols();

        gsMatrix<TT> paramBoundaryNodesFull;
        paramBoundaryNodesFull.setZero(patches.dim(),
paramBoundaryNodes.cols());
        gsMatrix<TT> paramFirstNodesFull = paramBoundaryNodesFull;

        int param = side.parameter();
        int index = 0;
        for (int s = 0; s < patches.dim(); s++)
        {
            if (s == side.direction())
            {
                if (param == 1) // east or north or back
                {
                    paramBoundaryNodesFull.middleRows(s, 1).setOnes();
                    TT nodeValue = basis->knot(s,
basis->knots(s).size() - basis->degree(s) - 2);

                    paramFirstNodesFull.middleRows(s,
1).setConstant(nodeValue);
                }
                else // param = 0 ... i.e. west or south or front
                {
                    TT nodeValue = basis->knot(s, basis->degree(s)
+ 1);
                    paramFirstNodesFull.middleRows(s,
1).setConstant(nodeValue);
                }

            }
            else
            {
                paramBoundaryNodesFull.middleRows(s, 1) =
paramBoundaryNodes.middleRows(index, 1);
                paramFirstNodesFull.middleRows(s, 1) =
paramBoundaryNodes.middleRows(index, 1);
                index++;
            }
        }

        gsMatrix<TT> physBoundaryNodesFull;
        gsMatrix<TT> physFirstNodesFull;
patches.patch(patchIndex).eval_into(paramBoundaryNodesFull,
physBoundaryNodesFull); // physical coordinates of paramBoundaryNodes
        patches.patch(patchIndex).eval_into(paramFirstNodesFull,
physFirstNodesFull); // physical coordinates of paramBoundaryNodes

        gsVector<TT> yPlusAtSide;
        yPlusAtSide.setZero(numOfPoints);
        for (int nodeIndex = 0; nodeIndex < numOfPoints; nodeIndex++)
        {
            TT distanceToWall = (physFirstNodesFull.col(nodeIndex) -
physBoundaryNodesFull.col(nodeIndex)).norm();
            yPlusAtSide[nodeIndex] = distanceToWall * u_tau / viscosity;
        }
        minYPlusOverSides[numSide] = yPlusAtSide.minCoeff();
    }
    return minYPlusOverSides.minCoeff() * viscosity / u_tau;
}
*/


