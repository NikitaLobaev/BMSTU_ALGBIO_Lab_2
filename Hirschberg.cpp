#include <algorithm>
#include <map>
#include "Matrix.h"
#include "NeedlemanWunsch.cpp"

namespace Lobaev::Hirschberg {
    
    template <class T, class L>
    std::pair<std::vector<T>, L> hirschberg_rec(const std::map<T, size_t> &matrix_map,
                                                const Lobaev::Math::Matrix<L> &matrix,
                                                L gap,
                                                const std::vector<T> &seq1,
                                                const std::vector<T> &seq2
    );

    template <class T, class L>
    std::pair<std::vector<T>, L> hirschberg(const std::map<T, size_t> &matrix_map,
                                            const Lobaev::Math::Matrix<L> &matrix,
                                            const L gap,
                                            const std::vector<T> &seq1,
                                            const std::vector<T> &seq2
    ) {
        if (!matrix.is_square()) {
            throw "hirschberg: only square similar matrix is allowed.";
        }
        if (matrix.rows_count() != matrix_map.size()) {
            throw "hirschberg: map of similar matrix is invalid, its size should be " +
                  std::to_string(matrix_map.size()) + ".";
        }

        return hirschberg_rec(matrix_map, matrix, gap, seq1, seq2);
    }
    
    template <class T, class L>
    std::vector<L> hirschberg_dp(const std::map<T, size_t> &matrix_map,
                                 const Lobaev::Math::Matrix<L> &matrix,
                                 L gap,
                                 const std::vector<T> &seq1,
                                 const std::vector<T> &seq2
    );
    
    template <class T, class L>
    std::pair<std::vector<T>, L> needleman_wunsch_unoptimized(const std::map<T, size_t> &matrix_map,
                                                              const Lobaev::Math::Matrix<L> &matrix,
                                                              L gap,
                                                              const std::vector<T> &seq1,
                                                              const std::vector<T> &seq2
    );
    
    template <class T, class L>
    std::pair<std::vector<T>, L> hirschberg_rec(const std::map<T, size_t> &matrix_map,
                                                const Lobaev::Math::Matrix<L> &matrix,
                                                const L gap,
                                                const std::vector<T> &seq1,
                                                const std::vector<T> &seq2
    ) {
        if (seq1.size() <= 1 || seq2.size() <= 1) { //допустимо решить такую маленькую задачу алгоритмом Нидлмана-Вунша, на асимтотическую сложность не повлияет
            return needleman_wunsch_unoptimized(matrix_map, matrix, gap, seq1, seq2);
        }
        
        std::vector<L> dp_left, dp_right;
        { //вычисление левого столбца матрицы динамического программирования (не использовал алгоритм Нидлмана-Вунша, т.к. он не возвращает этот столбец, переделывать его под это крайне неудобно)
            std::vector<T> seq1_(seq1.cbegin(), seq1.cbegin() + seq1.size() / 2);
            const std::vector<T> &seq2_ = seq2;
            dp_left = hirschberg_dp(matrix_map, matrix, gap, seq1_, seq2_);
        }
        { // -/- правого столбца -/-
            std::vector<T> seq1_(seq1.crbegin(), seq1.crbegin() + seq1.size() / 2);
            std::vector<T> seq2_(seq2.crbegin(), seq2.crend());
            dp_right = hirschberg_dp(matrix_map, matrix, gap, seq1_, seq2_);
            std::reverse(dp_right.begin(), dp_right.end());
        }
        
        if (dp_left.size() != dp_right.size() || dp_left.empty()) {
            throw "hirschberg: combine error.";
        }
        
        L seq2_max = dp_left[1] + dp_right[1];
        size_t seq2_max_index = 0;
        for (size_t i = 2; i < dp_left.size(); i++) { //поиск индекса, в котором сумма элементов левого и правого столбцов из матриц динамического программирования наименьшая
            if (seq2_max < dp_left[i] + dp_right[i]) {
                seq2_max = dp_left[i] + dp_right[i];
                seq2_max_index = i - 1;
            }
        }
        
        std::pair<std::pair<std::vector<T>, L>, std::pair<std::vector<T>, L>> rec_result;
        { //вычисление результата для левого верхнего квадранта
            const std::vector<T> seq1_(seq1.cbegin(), seq1.cbegin() + seq1.size() / 2);
            const std::vector<T> seq2_(seq2.cbegin(), seq2.cbegin() + seq2_max_index + 1);
            rec_result.first = hirschberg_rec(matrix_map, matrix, gap, seq1_, seq2_);
        }
        { // -/- для правого нижнего -/-
            const std::vector<T> seq1_(seq1.cbegin() + seq1.size() / 2, seq1.cend());
            const std::vector<T> seq2_(seq2.cbegin() + seq2_max_index + 1, seq2.cend());
            rec_result.second = hirschberg_rec(matrix_map, matrix, gap, seq1_, seq2_);
        }
        rec_result.first.first.insert(rec_result.first.first.end(),
                                      rec_result.second.first.cbegin(), rec_result.second.first.cend()); //результирующая строка будет храниться в rec_result.first.first
        return std::make_pair(rec_result.first.first,
                              rec_result.first.second + rec_result.second.second);
    }
    
    template <class T, class L>
    std::vector<L> hirschberg_dp(const std::map<T, size_t> &matrix_map,
                                 const Lobaev::Math::Matrix<L> &matrix,
                                 const L gap,
                                 const std::vector<T> &seq1,
                                 const std::vector<T> &seq2
    ) { //аналог алгоритма Нидлмана-Вунша; возвращает крайний столбец матрицы динамического программирования
        std::vector<L> dp(seq2.size() + 1);
        for (size_t i = 1; i < dp.size(); i++) {
            dp[i] = dp[i - 1] + gap;
        }
    
        for (size_t seq1_index = 0; seq1_index < seq1.size(); seq1_index++) {
            const size_t map_index_i = matrix_map.at(seq1[seq1_index]);
            
            L dp_prev = dp[0]; //для оптимизации реализованы такие нетривиальные записи в дополнительные переменные, чтобы хранить в памяти всего один std::vector<L> dp
            dp[0] += gap;
            for (size_t seq2_index = 0; seq2_index < seq2.size(); seq2_index++) {
                const size_t map_index_j = matrix_map.at(seq2[seq2_index]);
                
                L dp_new = dp_prev + matrix(map_index_i, map_index_j);
                dp_new = std::max(dp_new, dp[seq2_index] + gap);
                dp_new = std::max(dp_new, dp[seq2_index + 1] + gap);
                
                dp_prev = dp[seq2_index + 1];
                dp[seq2_index + 1] = dp_new;
            }
        }
    
        return dp;
    }
    
    template <class T, class L>
    std::pair<std::vector<T>, L> needleman_wunsch_unoptimized(const std::map<T, size_t> &matrix_map,
                                                              const Lobaev::Math::Matrix<L> &matrix,
                                                              const L gap,
                                                              const std::vector<T> &seq1,
                                                              const std::vector<T> &seq2
    ) { //аналог алгоритма Нидлмана-Вунша; возвращает крайний столбец матрицы динамического программирования
        std::vector<std::vector<L>> dp(seq1.size() + 1, std::vector<L>(seq2.size() + 1));
        for (size_t i = 1; i < dp.size(); i++) {
            dp[0][i] = dp[0][i - 1] + gap;
            dp[i][0] = dp[0][i];
        }
    
        for (size_t i = 1; i <= seq1.size(); i++) {
            for (size_t j = 1; j <= seq2.size(); j++) {
                const size_t map_index_i = matrix_map.at(seq1[i - 1]);
                const size_t map_index_j = matrix_map.at(seq2[j - 1]);
    
                dp[i][j] = dp[i - 1][j - 1] + matrix(map_index_i, map_index_j);
                dp[i][j] = std::max(dp[i][j], dp[i - 1][j] + gap);
                dp[i][j] = std::max(dp[i][j], dp[i][j - 1] + gap);
            }
        }
    
        L max = dp.back()[0];
        size_t max_index = 0;
        for (size_t i = 1; i < dp.back().size(); i++) { //поиск индекса, в котором сумма элементов левого и правого столбцов из матриц динамического программирования наименьшая
            if (max < dp.back()[i]) {
                max = dp.back()[i];
                max_index = i;
            }
        }
        
        std::pair<std::vector<T>, L> result;
        result.second = dp.back()[max_index];
        
        size_t i = max_index, j = dp.size() - 1;
        while (i > 0 && j > 0) {
            if (dp[i][j] == dp[i - 1][j] + gap) {
                i--;
            } else if (dp[i][j] == dp[i][j - 1] + gap) {
                j--;
            } else {
                i--;
                j--;
                result.first.emplace_back(seq1[i]);
            }
        }
    
        std::reverse(result.first.begin(), result.first.end());
    
        return result;
    }

}
