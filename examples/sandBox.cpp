#include <gismo.h>
#include <bitset>

using namespace gismo;

int main(int argc, char *argv[])
{
    int wanted = 4;

    gsMatrix<int> mat(3,2);
    mat << 1, 2, 3, 4, 5, 6;

    bool b = (mat.col(1).array() == wanted).any();
    std::string str = (b == 1) ? "true" : "false";

    gsInfo << "Is any element equal to " << wanted << "? - " << str << "\n";

    // ----------------------------------------------------------------------------

    // std::vector<gsMatrix<> > vec(3);
    // gsMatrix<> mat(2, 2);
    // mat << 1, 2, 3, 4;
    // vec[1] = mat;
    
    // gsInfo << "vec[0] = " << vec[0] << "\n";
    // gsInfo << "vec[1] = " << vec[1] << "\n";
    // gsInfo << "vec[2] = " << vec[2] << "\n";

    // ----------------------------------------------------------------------------

    // unsigned flag = NEED_VALUE;

    // gsInfo << "flag = " << flag << "\n";
    // gsInfo << "binary flag = " << std::bitset<10>(flag) << "\n";

    // flag = flag|NEED_MEASURE;

    // gsInfo << "flag = " << flag << "\n";
    // gsInfo << "binary flag = " << std::bitset<10>(flag) << "\n";

    // flag = flag|NEED_MEASURE;

    // gsInfo << "flag = " << flag << "\n";
    // gsInfo << "binary flag = " << std::bitset<10>(flag) << "\n";

    // unsigned flag1 = NEED_DERIV;

    // gsInfo << "flag & flag1 = " << (flag&flag1) << "\n";
    // gsInfo << "binary flag & flag1 = " << std::bitset<10>(flag&flag1) << "\n";

    // ----------------------------------------------------------------------------

    // gsKnotVector<> kv(0, 1, 3, 3);
    // gsTensorBSplineBasis<2> basis(kv, kv);

    // gsVector<real_t> point(2);
    // point << 0, 0; 

    // index_t basisID = 8;
    // gsMatrix<> supp = basis.support(basisID);
    // gsMatrix<index_t> elem = basis.elementSupport(basisID);
    // //size_t elemID = basis.elementIndex(point); // pada

    // gsInfo << "\nkv = " << kv << "\n";
    // gsInfo << "\nsupport " << basisID  << ":\n" << supp << "\n";
    // gsInfo << "\nelementSupport " << basisID  << ":\n" << elem << "\n";
    // //gsInfo << "elementIndex( " << point  << " ):\n" << elemID << "\n";

    // typename gsBasis<>::domainIter domIt = basis.makeDomainIterator(boundary::none);

    // while(domIt->good())
    // {
    //     bool inSupport = true; 

    //     for (index_t d = 0; d < 2; d++)
    //     {
    //         if ( (domIt->lowerCorner()[d] < supp(d,0)) ||  (domIt->upperCorner()[d] > supp(d,1)))
    //         {
    //             inSupport = false;
    //             break;
    //         }
    //     }

    //     gsInfo << "\nElement " << domIt->id() << " is in support of basis fcn " << basisID << ": " << inSupport;

    //     domIt->next();
    // }

    // gsInfo << "\n\n" << supp.middleCols(1,1);

}