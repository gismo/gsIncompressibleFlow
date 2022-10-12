/*
    Author(s): J. Egermaier
*/

#include <iostream>
#include <ctime>
#include <cstdlib>

#include <gismo.h>

using namespace gismo;

template<class T>
class uwbMG
{
public:
    uwbMG(int levels, gsOptionList optL)
    {
        m_levels = levels;
        m_fine = m_levels - 1;
        m_A.resize(m_levels); //system matrices
        m_M.resize(m_levels); //smoother matrices
        m_R.resize(m_fine); //transfer matrices
        m_r.resize(m_levels); //residual vector
        m_x.resize(m_levels); //solution (correction) vector
        m_basisU.resize(m_levels);
        m_basisP.resize(m_levels);
        m_opt = optL;
//        m_Rlin.resize(m_levels);
    }

    uwbMG(int levels, std::vector<gsSparseMatrix<T, RowMajor>> TM)
    {
        m_levels = levels;
        m_fine = m_levels - 1;
        m_A.resize(m_levels);
        m_M.resize(m_levels);
        m_R.resize(m_levels - 1);
        m_r.resize(m_levels);
        m_x.resize(m_levels);
        m_basisU.resize(m_levels);
        m_basisP.resize(m_levels);

        for (int i = 0; i < m_fine; i++){
            m_R[i] = TM[i];
            m_r[i].resize(m_R[i].cols());
            m_x[i].resize(m_R[i].cols());
            //gsInfo << "TM: " << TM[i].rows() << "x" << TM[i].cols() << "\n";
            //gsInfo << "m_R: " << m_R[i].rows() << "x" << m_R[i].cols() << "\n";
        }
    }

public:
    int m_levels, m_fine, usedIter;
    real_t m_gamma, relError;
    std::vector<gsSparseMatrix<T>> m_A, m_M, m_MT;
    std::vector<gsSparseMatrix<T, RowMajor>> m_R, m_RT;
    std::vector<gsVector<T>> m_r, m_x;
    std::vector<gsMultiBasis<T>> m_basisU, m_basisP;
    gsMatrix<T> m_rhs;
    gsMultiGridOp<>::Ptr m_mg;
    gsPreconditionerOp<>::Ptr m_smoother;
    gsOptionList m_opt;
    std::vector<gsLinearOperator<T>> m_Rlin;
    gsBoundaryConditions<T> m_BCU, m_BCP;
    std::vector<std::vector<gsMultiBasis<T>>> m_basisLevels;

public:

    void set_gamma(T gamma){ m_gamma = gamma; }
    void set_initialSolution(gsVector<T> x_input){ m_x[m_fine] = x_input; }
    void set_rhs(gsMatrix<T> rhs){ m_rhs = rhs; }
    void set_initialResiduum(gsVector<T> ir){ m_r[m_fine] = ir; }
    void set_smoothers(std::vector<gsSparseMatrix<T>> M){ m_M = M; }
    void set_optionList(gsOptionList OL){ m_opt = OL; }
    gsMultiBasis<T> get_fineBasisU(){ return m_basisU[m_fine]; }
    gsMultiBasis<T> get_fineBasisP(){ return m_basisP[m_fine]; }
    gsBoundaryConditions<T> get_BCU(){ return m_BCU; }
    gsBoundaryConditions<T> get_BCP(){ return m_BCP; }
    gsVector<T> get_correction(){ return m_x[m_fine]; }
    gsSparseMatrix<T> get_fineMatrix(){ return m_A[m_fine]; }
    gsMatrix<T> get_rhs(){ return m_rhs; }
    std::vector<gsSparseMatrix<T>> get_sysMatrices(){ return m_A; }
    std::vector<gsSparseMatrix<T, RowMajor>> get_transferMatrices(){ return m_R; }

    /*    void createMatrices(gsSparseMatrix<T> A)
    {
        m_A[m_fine] = A;
        for (int i = 0; i < m_fine; i++){
            int thisLevel = m_levels-2-i;
            int thisSize = m_R[thisLevel].cols();
            m_A[thisLevel].resize(thisSize, thisSize);
            //gsInfo << "rozmery : " << m_A[thisLevel].rows() << "x" << m_A[thisLevel].cols() << " = " << m_R[thisLevel].transpose().rows() << "x" << m_R[thisLevel].transpose().cols() << " * " << m_A[thisLevel+1].rows() << "x" <<m_A[thisLevel+1].cols() << " , " << m_R[thisLevel].rows() << "x" << m_R[thisLevel].cols() << "\n";
            m_A[thisLevel] = ((m_R[thisLevel].transpose()*m_A[thisLevel+1])*m_R[thisLevel]).pruned();
        }
    }
*/

    void createMatrices(gsSparseMatrix<T> A)
    {
        m_A[m_fine] = A;
        m_mg = gsMultiGridOp<>::make( m_A[m_fine], m_R );
        for (int i = 0; i < m_fine; i++){
            m_A[i] = m_mg->matrix(i);
        }
    }

    void createMatrices()
    {
        gsInfo << m_A.size() << m_R.size() << "\n";
        for (int i = 0; i < m_R.size(); i++){
            gsInfo << m_A[m_fine].rows() << "x" << m_A[m_fine].cols() << ", " << m_R[i].rows() << "x" << m_R[i].cols() << "\n";
        }
        m_mg = gsMultiGridOp<>::make( m_A[m_fine], m_R );
        for (int i = 0; i < m_fine; i++){
            m_A[i] = m_mg->matrix(i);
        }
    }

    void loadSystemAndBasis(std::string matrixName, gsVector<std::string> eraseName, bool saveToDat)
    {
        gsReadFile<real_t> (m_opt.getString("MG.matPath") + matrixName, m_A[m_fine]);
        //gsReadFile<real_t> (m_opt.getString("MG.matPath") + "step2D_visc_0-01_ref_4_0_0_deg_2_elev_0_st_blockA.xml", m_rhs);
        gsReadFile<real_t> (m_opt.getString("MG.matPath") + matrixName, m_rhs);

        util::string_replace(matrixName,".xml","");
        if (saveToDat){
            write_sparse_direct_to_dat(m_opt.getString("MG.matPath") + matrixName + ".dat" , m_A[m_fine]);
            gsFileData<real_t> fd;
            fd << m_rhs;
            fd.save(m_opt.getString("MG.matPath") + matrixName + "_rhs.xml");
        }

        for (int i = 0; i < eraseName.size(); i++){
            util::string_replace(matrixName,eraseName(i),"");
        }

        gsReadFile<real_t> (m_opt.getString("MG.domainPath") + matrixName + "discreteBasis_U.xml", m_basisU[m_fine]);
        gsReadFile<real_t> (m_opt.getString("MG.domainPath") + matrixName + "discreteBasis_P.xml", m_basisP[m_fine]);
        gsReadFile<real_t> (m_opt.getString("MG.domainPath") + matrixName + "BCU.xml", m_BCU);
        gsReadFile<real_t> (m_opt.getString("MG.domainPath") + matrixName + "BCP.xml", m_BCP);

    }

    void solve()
    {
        Vcycle(m_fine);
    }

    void setSmoother(int iL, std::string smoother)
    {
        if ( smoother == "Richardson" || smoother == "Ri" )
            m_smoother = makeRichardsonOp(m_mg->matrix(iL));
        else if ( smoother == "Jacobi" || smoother == "J" )
            m_smoother = makeJacobiOp(m_mg->matrix(iL));
        else if ( smoother == "GaussSeidel" || smoother == "GS" )
            m_smoother = makeGaussSeidelOp(m_mg->matrix(iL));
        else if ( smoother == "SubspaceCorrectedMassSmoother" || smoother == "SCMS" || smoother == "Hybrid" || smoother == "HYB" )
        {
            if (m_opt.getString("MG.BCtype") == "u"){
                if (m_basisLevels[0][iL].nBases() == 1){
                    m_smoother = gsPreconditionerFromOp<>::make(m_mg->underlyingOp(iL),gsPatchPreconditionersCreator<>::subspaceCorrectedMassSmootherOp(m_basisLevels[0][iL][0],m_BCU,m_opt.getGroup("MG"),m_opt.getReal("MG.scaling")));
                } else {
                    m_smoother = setupSubspaceCorrectedMassSmoother( m_mg->matrix(iL), m_basisLevels[0][iL], m_BCU, m_opt.getGroup("MG") );
                }
            } else {
                if (m_basisLevels[0][iL].nBases() == 1){
                    m_smoother = gsPreconditionerFromOp<>::make(m_mg->underlyingOp(iL),gsPatchPreconditionersCreator<>::subspaceCorrectedMassSmootherOp(m_basisLevels[0][iL][0],m_BCP,m_opt.getGroup("MG"),m_opt.getReal("MG.scaling")));
                } else {
                    m_smoother = setupSubspaceCorrectedMassSmoother( m_mg->matrix(iL), m_basisLevels[0][iL], m_BCP, m_opt.getGroup("MG") );
                }
            }

            if ( smoother == "Hybrid" || smoother == "HYB" )
            {
                m_smoother->setOptions( m_opt.getGroup("MG") );
                m_smoother = gsCompositePrecOp<>::make( makeGaussSeidelOp(m_mg->matrix(iL)), m_smoother );
            }
        } else {
            gsInfo << "\n\nThe chosen smoother is unknown.\n\nKnown are:\n  Richardson (Ri)\n  Jacobi (J)\n  GaussSeidel (GS)"
                      "\n  SubspaceCorrectedMassSmoother (SCMS)\n  Hybrid (HYB)\n\n";
        }
    }

    void solveAsPrecond(std::string smoother, std::string method, bool print = false)
    {
        m_mg->setOptions( m_opt.getGroup("MG") );
        for (index_t i = 1; i < m_mg->numLevels(); ++i){
            setSmoother(i,smoother);
            m_smoother->setOptions( m_opt.getGroup("MG") );
            // Handle the extra-smooth option. On the finest grid level, there is nothing to handle.
            if (m_opt.getSwitch("MG.Extrasmooth") && i < m_mg->numLevels()-1 ){
                m_smoother->setNumOfSweeps( 1 << (m_mg->numLevels()-1-i) );
                m_smoother = gsPreconditionerFromOp<>::make(m_mg->underlyingOp(i),m_smoother);
            }
            m_mg->setSmoother(i, m_smoother);
            //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            //m_Rlin[i] = *m_mg->underlyingOp(i);
            //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        }

        gsMatrix<> xMat, errorHistory;
        //xMat.setRandom( m_mg->matrix(m_fine).rows(), 1 );
        xMat.resize(m_x[m_fine].size(), 1);
        xMat.col(0) = m_x[m_fine];

        if (method=="CG")
            gsConjugateGradient<>( m_A[m_fine], m_mg )
                .setOptions( m_opt.getGroup("Solver") )
                .solveDetailed( m_rhs, xMat, errorHistory );
        else if (method=="GM")
            gsGradientMethod<>( m_A[m_fine], m_mg )
                .setOptions( m_opt.getGroup("Solver") )
                .solveDetailed( m_rhs, xMat, errorHistory );
        else {
            gsInfo << "\n\nThe chosen iterative solver is unknown.\n\nKnown are:\n  conjugate gradient (CG)\n  direct (d)\n\n";
        }

        const index_t iter = errorHistory.rows()-1;
        if(print){
            const bool success = errorHistory(iter,0) < m_opt.getReal("Solver.Tolerance");
            if (success)
                gsInfo << "Reached desired tolerance after " << iter << " iterations:\n";
            else
                gsInfo << "Did not reach desired tolerance after " << iter << " iterations:\n";

            if (errorHistory.rows() < 20)
                gsInfo << errorHistory.transpose() << "\n\n";
            else
                gsInfo << errorHistory.topRows(5).transpose() << " ... " << errorHistory.bottomRows(5).transpose()  << "\n\n";
        }
        m_x[m_fine] = xMat.col(0);
        usedIter = iter;
        relError = errorHistory(iter);
    }

    void set_smoothers()
    {
        for(int i = 0; i < m_levels; i++){
            gsSparseMatrix<T> L, U, D, Dinv;
            D = m_A[i].diagonal().asDiagonal();
            //U = m_A[i].triangularView<Upper>() - D;
            //L = m_A[i].triangularView<Lower>() - D;
            Dinv = (m_A[i].diagonal().asDiagonal()).inverse();
            m_M[i] = -Dinv*(L+U);
        }
    }

    void createTransferMatrices (std::vector<gsMultiBasis<T>> MB, std::vector<gsBoundaryConditions<T>> BC)
    {
        std::vector< std::vector< gsSparseMatrix<T, RowMajor> > > vecTRM;
        m_basisLevels.resize(MB.size());
        for (unsigned i = 0; i < MB.size(); i++){
            std::vector< gsSparseMatrix<T, RowMajor> > trM;
            gsGridHierarchy<>::buildByCoarsening(give(MB[i]), BC[i], m_opt.getGroup("MG"))
                .moveMultiBasesTo( m_basisLevels[i] )
                .moveTransferMatricesTo(trM)
                .clear();
            vecTRM.push_back(trM);
        }
        for (unsigned j = 0; j < vecTRM[0].size(); j++){
            std::vector<gsSparseMatrix<T, RowMajor>> vecM;
            for (unsigned i = 0; i < vecTRM.size(); i++){
                vecM.push_back(vecTRM[i][j]);
            }
            m_R[j] = combineMatrices(vecM);
        }
    }

    gsPreconditionerOp<>::Ptr setupSubspaceCorrectedMassSmoother(const gsSparseMatrix<>& matrix, const gsMultiBasis<>& mb, const gsBoundaryConditions<>& bc, const gsOptionList& opt)
    {
        const short_t dim = mb.topology().dim();

        // Setup dof mapper
        gsDofMapper dm;
        mb.getMapper((dirichlet::strategy)opt.askInt("DirichletStrategy",11),(iFace::strategy)opt.askInt("InterfaceStrategy", 1),bc,dm,0);
        const index_t nTotalDofs = dm.freeSize();

        // Decompose the whole domain into components
        std::vector< std::vector<patchComponent> > components = mb.topology().allComponents(true);
        const index_t nr_components = components.size();

        // Setup Dirichlet boundary conditions
        gsBoundaryConditions<> dir_bc;
        for( index_t ps=0; ps < 2*dim; ++ps )
            dir_bc.addCondition( 0, 1+ps, condition_type::dirichlet, NULL );

        // Setup transfer matrices and local preconditioners
        std::vector< gsSparseMatrix<real_t,RowMajor> > transfers;
        transfers.reserve(nr_components);

        std::vector< gsLinearOperator<>::Ptr > ops;
        ops.reserve(nr_components);

        for (index_t i=0; i<nr_components; ++i)
        {
            gsMatrix<unsigned> indices;
            std::vector<gsBasis<>::uPtr> bases = mb.componentBasis_withIndices(components[i],dm,indices,true);

            index_t sz = indices.rows();
            gsSparseEntries<> se;
            se.reserve(sz);
            for (index_t j=0; j<sz; ++j)
                se.add(indices(j,0),j,real_t(1));
            gsSparseMatrix<real_t,RowMajor> transfer(nTotalDofs,sz);
            transfer.setFrom(se);

            if (sz>0)
            {
                if (bases[0]->dim() == dim)
                {
                    GISMO_ASSERT ( bases.size() == 1, "Only one basis is expected for each patch." );
                    ops.push_back(gsPatchPreconditionersCreator<>::subspaceCorrectedMassSmootherOp(*(bases[0]),dir_bc,gsOptionList(),opt.getReal("Scaling")));
                }
                else
                {
                    gsSparseMatrix<> mat = transfer.transpose() * matrix * transfer;
                    ops.push_back( makeSparseCholeskySolver(mat) );
                }
                transfers.push_back(give(transfer));
            }
        }

        return gsPreconditionerFromOp<>::make(makeMatrixOp(matrix),gsAdditiveOp<>::make(transfers, ops));
    }

    void Vcycle(int j)
    {
        if (j > 0){
            m_x[j].resize(m_A[j].rows());
            m_x[j].setZero(m_A[j].rows());

            m_x[j] += m_gamma*m_M[j]*m_r[j];
            m_r[j] -= m_A[j]*m_x[j];
            m_r[j-1] = m_R[j-1].transpose()*m_r[j];
            Vcycle(j-1);
            m_x[j] += m_R[j-1]*m_x[j-1];
            m_r[j] -= m_A[j]*m_x[j];

            m_x[j] += m_gamma*m_M[j].transpose()*m_r[j];
        } else {
            typename gsSparseSolver<T>::LU directSolver;
            directSolver.compute(m_A[j]);
            m_x[j] = directSolver.solve(m_r[j]);
        }
    }

    void addAsBlock(const gsSparseMatrix<T>& mat, index_t row, index_t col, gsSparseMatrix<T>& result)
    {
        GISMO_ASSERT( row + mat.rows() <= result.rows(), "Dimensions to not agree" );
        GISMO_ASSERT( col + mat.cols() <= result.cols(), "Dimensions to not agree" );
        //TODO: this is not really good in terms of efficiency
        for (index_t k=0; k < mat.outerSize(); ++k)
            for (typename gsSparseMatrix<T>::InnerIterator it(mat,k); it; ++it)
                result(row+it.row(),col+it.col()) += it.value();
    }

    void addAsBlock(const gsSparseMatrix<T,RowMajor>& mat, index_t row, index_t col, gsSparseMatrix<T,RowMajor>& result)
    {
        GISMO_ASSERT( row + mat.rows() <= result.rows(), "Dimensions to not agree" );
        GISMO_ASSERT( col + mat.cols() <= result.cols(), "Dimensions to not agree" );
        //TODO: this is not really good in terms of efficiency
        for (index_t k=0; k < mat.outerSize(); ++k)
            for (typename gsSparseMatrix<T,RowMajor>::InnerIterator it(mat,k); it; ++it)
                result(row+it.row(),col+it.col()) += it.value();
    }

    gsSparseMatrix<T, RowMajor> combineMatrices(std::vector<gsSparseMatrix<T, RowMajor>> vecMat)
    {
        gsSparseMatrix<T, RowMajor> output;
        int rows = 0;
        int cols = 0;
        for (unsigned i = 0; i < vecMat.size(); i++){
            rows += vecMat[i].rows();
            cols += vecMat[i].cols();
        }
        output.resize(rows,cols);

        int row_cnt = 0;
        int col_cnt = 0;
        for (unsigned i = 0; i < vecMat.size(); i++){
            //gsInfo << "rows = " << rows << ", cols = " << cols << ", row_cnt = " << row_cnt << ", col_cnt = " << col_cnt << ", i = " << i << ", vecMat[i].rows() = " << vecMat[i].rows() << ", vecMat[i].cols() = " << vecMat[i].cols() << "\n";
            addAsBlock(vecMat[i], row_cnt, col_cnt, output);
            row_cnt = vecMat[i].rows();
            col_cnt = vecMat[i].cols();
        }
        return output;
    }

    gsSparseMatrix<T> combineMatrices(std::vector<gsSparseMatrix<T>> vecMat)
    {
        gsSparseMatrix<T> output;
        int rows = 0;
        int cols = 0;
        for (unsigned i = 0; i < vecMat.size(); i++){
            rows += vecMat[i].rows();
            cols += vecMat[i].cols();
        }
        output.resize(rows,cols);

        int row_cnt = 0;
        int col_cnt = 0;
        for (unsigned i = 0; i < vecMat.size(); i++){
            //gsInfo << "rows = " << rows << ", cols = " << cols << ", row_cnt = " << row_cnt << ", col_cnt = " << col_cnt << ", i = " << i << ", vecMat[i].rows() = " << vecMat[i].rows() << ", vecMat[i].cols() = " << vecMat[i].cols() << "\n";
            addAsBlock(vecMat[i], row_cnt, col_cnt, output);
            row_cnt = vecMat[i].rows();
            col_cnt = vecMat[i].cols();
        }
        return output;
    }
};
