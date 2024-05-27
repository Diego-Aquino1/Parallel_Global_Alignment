#include <iostream>
#include <vector>
#include <fstream>
#include <functional>
#include <omp.h>

class NeedlemanWunsch {
private:
    std::string seq1;
    std::string seq2;
    int match_score;
    int mismatch_score;
    int gap_penalty;
    int m;
    int n;
    std::vector<std::vector<int>> score_matrix;
    std::vector<std::vector<std::vector<std::pair<int, int>>>> traceback_matrix;
    std::vector<std::pair<std::string, std::string>> alignments;

public:
    NeedlemanWunsch(std::string seq1, std::string seq2, int match_score = 1, int mismatch_score = -1, int gap_penalty = -2) :
        seq1(seq1), seq2(seq2), match_score(match_score), mismatch_score(mismatch_score), gap_penalty(gap_penalty) {
        m = seq1.length();
        n = seq2.length();
        score_matrix.resize(m + 1, std::vector<int>(n + 1, 0));
        traceback_matrix.resize(m + 1, std::vector<std::vector<std::pair<int, int>>>(n + 1, std::vector<std::pair<int, int>>()));
    }

    void align() {
        double start_time = omp_get_wtime();
        initialize_matrices();
        fill_matrices();
        //traceback(m, n, "", "");
        write_results();
        double end_time = omp_get_wtime();
        double elapsed_time = (end_time - start_time) * 1000.0;
        std::cout << "Tiempo de ejecución: " << elapsed_time << " milisegundos\n";
    }

private:

    std::string removeSpaces(std::string str) {

        str.erase(remove(str.begin(), str.end(), ' '), str.end());
        return str;

    }

    void initialize_matrices() {
        for (int i = 1; i <= m; ++i) {
            score_matrix[i][0] = score_matrix[i - 1][0] + gap_penalty;
            traceback_matrix[i][0].push_back(std::make_pair(i - 1, 0));
        }
        for (int j = 1; j <= n; ++j) {
            score_matrix[0][j] = score_matrix[0][j - 1] + gap_penalty;
            traceback_matrix[0][j].push_back(std::make_pair(0, j - 1));
        }
    }

    void fill_matrices() {
        std::vector<std::vector<std::vector<std::pair<int, int>>>> local_traceback_matrix(m + 1, std::vector<std::vector<std::pair<int, int>>>(n + 1));

        #pragma omp parallel for collapse(2)
        for (int i = 1; i <= m; ++i) {
            for (int j = 1; j <= n; ++j) {
                int match = score_matrix[i - 1][j - 1] + (seq1[i - 1] == seq2[j - 1] ? match_score : mismatch_score);
                int delete_op = score_matrix[i - 1][j] + gap_penalty;
                int insert_op = score_matrix[i][j - 1] + gap_penalty;
                int max_score = std::max(std::max(match, delete_op), insert_op);
                score_matrix[i][j] = max_score;

                if (match == max_score) {
                    local_traceback_matrix[i][j].push_back(std::make_pair(i - 1, j - 1));
                }
                if (delete_op == max_score) {
                    local_traceback_matrix[i][j].push_back(std::make_pair(i - 1, j));
                }
                if (insert_op == max_score) {
                    local_traceback_matrix[i][j].push_back(std::make_pair(i, j - 1));
                }
            }
        }

        // Merge local_traceback_matrix into traceback_matrix
        #pragma omp parallel for collapse(2)
        for (int i = 1; i <= m; ++i) {
            for (int j = 1; j <= n; ++j) {
                traceback_matrix[i][j].insert(traceback_matrix[i][j].end(),
                                            local_traceback_matrix[i][j].begin(),
                                            local_traceback_matrix[i][j].end());
            }
        }
    }

    void traceback(int i, int j, std::string alignment1, std::string alignment2) {
        if (i == 0 && j == 0) {
            alignments.push_back(std::make_pair(std::string(alignment1.rbegin(), alignment1.rend()), std::string(alignment2.rbegin(), alignment2.rend())));
            return;
        }
        for (auto& prev : traceback_matrix[i][j]) {
            int prev_i = prev.first;
            int prev_j = prev.second;
            if (prev == std::make_pair(i - 1, j - 1)) {
                traceback(prev_i, prev_j, alignment1 + seq1[i - 1], alignment2 + seq2[j - 1]);
            }
            else if (prev == std::make_pair(i - 1, j)) {
                traceback(prev_i, prev_j, alignment1 + seq1[i - 1], alignment2 + '-');
            }
            else if (prev == std::make_pair(i, j - 1)) {
                traceback(prev_i, prev_j, alignment1 + '-', alignment2 + seq2[j - 1]);
            }
        }
    }

    void write_results() {
        std::ofstream file("D:/UNSA/5/Bioinformatica/Parallel_Global_Alignment/parallel_alignment_results.txt");
        if (file.is_open()) {
            file << "Score Matrix:\n";
            for (auto& row : score_matrix) {
                for (auto& val : row) {
                    file << val << "\t";
                }
                file << "\n";
            }
            int final_score = score_matrix[m][n];
            file << "\nFinal Score: " << final_score << "\n";

            int num_alignments = alignments.size();
            file << "Number of Optimal Alignments: " << num_alignments << "\n";

            file << "\nOptimal Alignments:\n";
            for (size_t i = 0; i < alignments.size(); ++i) {
                file << "Alignment " << i + 1 << ":\n";
                file << alignments[i].first << "\n";
                file << alignments[i].second << "\n\n";
            }
            file.close();
            std::cout << "Results written to parallel_alignment_results.txt\n";
        }
        else {
            std::cerr << "Unable to open file for writing.\n";
        }
    }
};

int main() {
    omp_set_num_threads(12);
    // std::string seq1 = "ACTGATTCA";
    // std::string seq2 = "ACGCATCA";

    std::string seq1 = "atggaagcaatatcactgatgactatactactggtggtaacaacaagtaatgcagacaaaatctgcatcggtcaccaatcaacaaattccacggaaactgtagacacgctaacagaaacaaatgttcctgtaacacaagccaaagaattgctccacacagaacacaatgggatgctatgtgcaacaaatctgggacgtcctcttatcctagacacatgcaccattgaaggactgatctatggcaacccatcttgtgacatgctgttaggaggaagggaatggtcctacatcgtcgaaagaccgtccgcagtaaatggaacatgctaccctggaaatgtagaaaacctagaggaacttagaacactttttagctcctctagttcttaccaaagagtccaactctttccagactcaatctggaatgtgacttacactgggacaagcaaatcatgttcagattcattctataggaatatgagatggttaactcaaaagaatgggggttatccaattcaagatgcccagtacacaaacaataggggaaaggacattcttttcgtgtggggcatacatcatccaccaaccgatactgcacagacgaatttatatacaaggaccgacacaacaacaagtgtaacaacggagactttagataggaccttcaaaccattgatagggccaaggccccttgtcaatggtctaattggaagaattaattactattggtcggtactaaaaccaggccaaacgttgcgagtgagatcaaatggaaatctaattgctccatggtttggacatgttctctcaggtgagagccatgtgagaatcctgagaactgatttaagcagcggtaattgtgtggtacaatgccagactgaaaaaggtggcctaaacagtacaatgccatttcacaacatcagcaaatatgcttttgggacctgtcccaaatatattggagtcaagagtctcaaactggcaattggccttagaaacgtacatgctaggtcaagtagaggactattcggagcgatagctggattcatagaaggaggttggccaggactagtcgccggttggtat";
    std::string seq2 = "attaaaggtttataccttcccaggtaacaaaccaaccaactttcgatctcttgtagatctgttctctaaacgaactttaaaatctgtgtggctgtcactcggctgcatgcttagtgcactcacgcagtataattaataactaattactgtcgttgacaggacacgagtaactcgtctatcttctgcaggctgcttacggtttcgtccgtgttgcagccgatcatcagcacatctaggtttcgtccgggtgtgaccgaaaggtaagatggagagccttgtccctggtttcaacgagaaaacacacgtccaactcagtttgcctgttttacaggttcgcgacgtgctcgtacgtggctttggagactccgtggaggaggtcttatcagaggcacgtcaacatcttaaagatggcacttgtggcttagtagaagttgaaaaaggcgttttgcctcaacttgaacagccctatgtgttcatcaaacgttcggatgctcgaactgcacctcatggtcatgttatggttgagctggtagcagaactcgaaggcattcagtacggtcgtagtggtgagacacttggtgtccttgtccctcatgtgggcgaaataccagtggcttaccgcaaggttcttcttcgtaagaacggtaataaaggagctggtggccatagttacggcgccgatctaaagtcatttgacttaggcgacgagcttggcactgatccttatgaagattttcaagaaaactggaacactaaacatagcagtggtgttacccgtgaactcatgcgtgagcttaacggaggggcatacactcgctatgtcgataacaacttctgtggccctgatggctaccctcttgagtgcattaaagaccttctagcacgtgctggtaaagcttcatgcactttgtccgaacaactggactttattgacactaagaggggtgtatactgctgccgtgaacatgagcatgaaattgcttggtacacggaacgttctgaaaagagctatgaattgcagacaccttttgaaattaaattggcaaagaaatttgacaccttcaatggggaatgtccaaa";

    std::cout << "COMPILO";

    NeedlemanWunsch nw(seq1, seq2);
    nw.align();
    return 0;
}
