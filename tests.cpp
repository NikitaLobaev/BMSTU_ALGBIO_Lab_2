#include "gtest/gtest.h"
#include "Hirschberg.cpp"

using namespace Lobaev::Math;
using namespace Lobaev::Hirschberg;

const Matrix<long> matrix_blosum62({
                                   { 4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0, -2, -1,  0, -4},
                                   {-1,  5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1, -3, -2, -1, -1, -3, -2, -3, -1,  0, -1, -4},
                                   {-2,  0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3,  3,  0, -1, -4},
                                   {-2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3,  4,  1, -1, -4},
                                   { 0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1, -3, -3, -2, -4},
                                   {-1,  1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2,  0,  3, -1, -4},
                                   {-1,  0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -3, -2, -2,  1,  4, -1, -4},
                                   { 0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3, -1, -2, -1, -4},
                                   {-2,  0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3,  0,  0, -1, -4},
                                   {-1, -3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,  0, -3, -2, -1, -3, -1,  3, -3, -3, -1, -4},
                                   {-1, -2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,  0, -3, -2, -1, -2, -1,  1, -4, -3, -1, -4},
                                   {-1,  2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1, -3, -1,  0, -1, -3, -2, -2,  0,  1, -1, -4},
                                   {-1, -1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,  0, -2, -1, -1, -1, -1,  1, -3, -1, -1, -4},
                                   {-2, -3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,  6, -4, -2, -2,  1,  3, -1, -3, -3, -1, -4},
                                   {-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1, -4, -3, -2, -2, -1, -2, -4},
                                   { 1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1, -2, -1,  4,  1, -3, -2, -2,  0,  0,  0, -4},
                                   { 0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0, -1, -1,  0, -4},
                                   {-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,  1, -4, -3, -2, 11,  2, -3, -4, -3, -2, -4},
                                   {-2, -2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,  3, -3, -2, -2,  2,  7, -1, -3, -2, -1, -4},
                                   { 0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3, -1,  4, -3, -2, -1, -4},
                                   {-2, -1,  3,  4, -3,  0,  1, -1,  0, -3, -4,  0, -3, -3, -2,  0, -1, -4, -3, -3,  4,  1, -1, -4},
                                   {-1,  0,  0,  1, -3,  3,  4, -2,  0, -3, -3,  1, -1, -3, -1,  0, -1, -3, -2, -2,  1,  4, -1, -4},
                                   { 0, -1, -1, -1, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -2,  0,  0, -2, -1, -1, -1, -1, -1, -4},
                                   {-4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4,  1},
                                   });
const std::map<char, size_t> map_blosum62{
    {'A', 0},
    {'R', 1},
    {'N', 2},
    {'D', 3},
    {'C', 4},
    {'Q', 5},
    {'E', 6},
    {'G', 7},
    {'H', 8},
    {'I', 9},
    {'L', 10},
    {'K', 11},
    {'M', 12},
    {'F', 13},
    {'P', 14},
    {'S', 15},
    {'T', 16},
    {'W', 17},
    {'Y', 18},
    {'V', 19},
    {'B', 20},
    {'Z', 21},
    {'X', 22},
    {'*', 23},
};

TEST(hirschberg, test_1_ok) {
    const long gap = -10;
    const std::string seq1 = "MHSKVTIICIRFLFWFLLLCMLIGKSHTEDDIIIATKNGKVRGMNLTVFGGTVTAFLGIP"
                             "YAQPPLGRLRFKKPQSLTKWSDIWNATKYANSCCQNIDQSFPGFHGSEMWNPNTDLSEDC"
                             "LYLNVWIPAPKPKNATVLIWIYGGGFQTGTSSLHVYDGKFLARVERVIVVSMNYRVGALG"
                             "FLALPGNPEAPGNMGLFDQQLALQWVQKNIAAFGGNPKSVTLFGESAGAASVSLHLLSPG"
                             "SHSLFTRAILQSGSFNAPWAVTSLYEARNRTLNLAKLTGCSRENETEIIKCLRNKDPQEI"
                             "LLNEAFVVPYGTPLSVNFGPTVDGDFLTDMPDILLELGQFKKTQILVGVNKDEGTAFLVY"
                             "GAPGFSKDNNSIITRKEFQEGLKIFFPGVSEFGKESILFHYTDWVDDQRPENYREALGDV"
                             "VGDYNFICPALEFTKKFSEWGNNAFFYYFEHRSSKLPWPEWMGVMHGYEIEFVFGLPLER"
                             "RDNYTKAEEILSRSIVKRWANFAKYGNPNETQNNSTSWPVFKSTEQKYLTLNTESTRIMT"
                             "KLRAQQCRFWTSFFPKVLEMTGNIDEAEWEWKAGFHRWNNYMMDWKNQFNDYTSKKESCV"
                             "GL";
    const std::string seq2 = "MQTQHTKVTQTHFLLWILLLCMPFGKSHTEEDFIITTKTGRVRGLSMPVLGGTVTAFLGI"
                             "PYAQPPLGSLRFKKPQPLNKWPDIHNATQYANSCYQNIDQAFPGFQGSEMWNPNTNLSED"
                             "CLYLNVWIPVPKPKNATVMVWIYGGGFQTGTSSLPVYDGKFLARVERVIVVSMNYRVGAL"
                             "GFLAFPGNPDAPGNMGLFDQQLALQWVQRNIAAFGGNPKSITIFGESAGAASVSLHLLCP"
                             "QSYPLFTRAILESGSSNAPWAVKHPEEARNRTLTLAKFTGCSKENEMEMIKCLRSKDPQE"
                             "ILRNERFVLPSDSILSINFGPTVDGDFLTDMPHTLLQLGKVKKAQILVGVNKDEGTAFLV"
                             "YGAPGFSKDNDSLITRKEFQEGLNMYFPGVSRLGKEAVLFYYVDWLGEQSPEVYRDALDD"
                             "VIGDYNIICPALEFTKKFAELENNAFFYFFEHRSSKLPWPEWMGVMHGYEIEFVFGLPLG"
                             "RRVNYTRAEEIFSRSIMKTWANFAKYGHPNGTQGNSTMWPVFTSTEQKYLTLNTEKSKIY"
                             "SKLRAPQCQFWRLFFPKVLEMTGDIDETEQEWKAGFHRWSNYMMDWQNQFNDYTSKKESC"
                             "TAL";
    
    const std::pair<std::vector<char>, long> result = hirschberg(map_blosum62, matrix_blosum62, gap,
                                                                 std::vector<char>(seq1.cbegin(), seq1.cend()),
                                                                 std::vector<char>(seq2.cbegin(), seq2.cend()));
    
    const long expected_score = 2373;
    
    ASSERT_EQ(result.second, expected_score);
}
