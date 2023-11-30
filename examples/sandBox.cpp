#include <gismo.h>
#include <bitset>

using namespace gismo;

int main(int argc, char *argv[])
{
    std::vector<gsMatrix<> > vec(3);
    gsMatrix<> mat(2, 2);
    mat << 1, 2, 3, 4;
    vec[1] = mat;
    
    gsInfo << "vec[0] = " << vec[0] << "\n";
    gsInfo << "vec[1] = " << vec[1] << "\n";
    gsInfo << "vec[2] = " << vec[2] << "\n";


    // unsigned flag = NEED_VALUE;

    // gsInfo << "flag = " << flag << "\n";
    // gsInfo << "binary flag = " << std::bitset<10>(flag) << "\n";

    // flag = flag|NEED_MEASURE;

    // gsInfo << "flag = " << flag << "\n";
    // gsInfo << "binary flag = " << std::bitset<10>(flag) << "\n";

    // flag = flag|NEED_MEASURE;

    // gsInfo << "flag = " << flag << "\n";
    // gsInfo << "binary flag = " << std::bitset<10>(flag) << "\n";


    // gsKnotVector<> kv(0, 1, 3, 3);
    // gsTensorBSplineBasis<2> basis(kv, kv);

    // index_t basisID = 8;
    // gsMatrix<> supp = basis.support(basisID);
    // gsMatrix<index_t> elem = basis.elementSupport(basisID);

    // gsInfo << "kv = " << kv << "\n";
    // gsInfo << "support " << basisID  << ":\n" << supp << "\n\n";
    // gsInfo << "elementSupport " << basisID  << ":\n" << elem << "\n";

}