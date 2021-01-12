// fluxRelated.h

// Erin Sam Joe | NITK '23 | Mechanical Engg. Major | Electroncs & Electrical Engg. Minor


#include <vector>
#include <complex>
#include <algorithm>
#include <Dense>
#include <Eigenvalues>

// Prototypes
double eigvalMaxMinF(const std::vector<double>&, bool=false);
double eigvalMaxMinG(const std::vector<double> U, bool=false);
void fluxF(const std::vector<double>& U, std::vector<double>& flux);
void fluxG(const std::vector<double>& U, std::vector<double>& flux);
