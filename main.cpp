#include <iostream>
#include "Hirschberg.cpp"
#include "IO.cpp"

typedef std::pair<std::string, std::string> DNA; //<dna_id, dna_string>

void read_dnas(std::istream&, std::vector<DNA>&);

const std::map<char, size_t> default_matrix_map{
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

//TODO
const char *usage = "Usage: lab1 (-m | --matrix) <input matrix filename> [(-i | --input) <input filename>] [(-o | --output) <output filename>] [(-g | --gap) <gap>]";

int main(int argc, char **argv) {
    std::string input_filename, output_filename, matrix_filename;
    long gap = -2;

    if (argc % 2 == 0) {
        std::cerr << usage << std::endl;
        return 1;
    }

    for (size_t i = 1; i + 1 < argc; i += 2) {
        std::string arg(argv[i]);

        if (arg == "-i" || arg == "--input") {
            input_filename = argv[i + 1];
        } else if (arg == "-o" || arg == "--output") {
            output_filename = argv[i + 1];
        } else if (arg == "-g" || arg == "--gap") {
            gap = std::stol(argv[i + 1]);
        } else if (arg == "-m" || arg == "--matrix") {
            matrix_filename = argv[i + 1];
        } else {
            std::cerr << usage << std::endl;
            return 1;
        }
    }

    if (!input_filename.empty()) {
        auto file = std::freopen(input_filename.c_str(), "r", stdin);
        if (!file) {
            std::cerr << "Input file doesn't exist." << std::endl << usage << std::endl;
            return 1;
        }
    }

    if (!output_filename.empty()) {
        auto file = std::freopen(output_filename.c_str(), "w", stdout);
        if (!file) {
            std::cerr << "Output file doesn't exist." << std::endl << usage << std::endl;
            return 1;
        }
    }

    std::ifstream matrix_ifstream(matrix_filename);
    if (!matrix_ifstream) {
        std::cerr << "Input matrix file doesn't exist." << std::endl << usage << std::endl;
        return 1;
    }

    auto matrix = Lobaev::IO::read_matrix<long>(matrix_ifstream);

    matrix_ifstream.close();

    std::vector<DNA> dnas;
    read_dnas(std::cin, dnas);

    if (dnas.size() != 2) {
        std::cerr << "Only 2 dna's in input file are allowed." << std::endl << usage << std::endl;
        return 1;
    }

    std::pair<std::vector<char>, long> result = Lobaev::Hirschberg::hirschberg<char>(default_matrix_map, matrix, gap,
               std::vector<char>(dnas[0].second.cbegin(), dnas[0].second.cend()),
               std::vector<char>(dnas[1].second.cbegin(), dnas[1].second.cend()));
    
    for (const char &c : result.first) {
        std::cout << c;
    }
    std::cout << std::endl;
    std::cout << "Score: " << result.second << std::endl;

    return 0;
}

void read_dnas(std::istream &in, std::vector<DNA> &dnas) {
    std::string cur_line, cur_dna_buf;
    while (std::getline(in, cur_line)) {
        if (cur_line.empty()) {
            continue;
        }

        if (cur_line[0] == '>') {
            if (!cur_dna_buf.empty()) {
                dnas.back().second = cur_dna_buf;
                cur_dna_buf.clear();
            }

            const size_t id_start_index = cur_line.find('|');
            const size_t id_length = cur_line.substr(id_start_index + 1).find('|');
            const std::string id = cur_line.substr(id_start_index + 1, id_length);
            dnas.emplace_back(id, "");
        } else {
            cur_dna_buf += cur_line;
        }
    }

    if (!cur_dna_buf.empty()) {
        dnas.back().second = cur_dna_buf;
    }
}
