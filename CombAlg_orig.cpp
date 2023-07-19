#include <iostream>
#include <vector>
#include <utility>
#include <algorithm>
#include <fstream>
#include <ctime>
#include <set>
#include <boost/dynamic_bitset.hpp>

int max_level = -100;
int r = 10000000;

bool comp(std::pair<std::vector<boost::dynamic_bitset<>>, int> &pat1, std::pair<std::vector<boost::dynamic_bitset<>>, int> &pat2)
{
    return pat1.first.size() > pat2.first.size();
}

int find_object(std::vector<boost::dynamic_bitset<>> &X0, std::vector<boost::dynamic_bitset<>> &X1)
{
    int n = X0[0].size();
    int m0 = X0.size();
    int m1 = X1.size();
    std::vector<int> dist(m1);
    for (int i = 0; i < m1; i++) {
        int d = 0;
        for (int j = 0; j < n; j++) {
            int feat_d = 0;
            for (int k = 0; k < m0; k++) {
                if (X0[k][j] != X1[i][j]) {
                    feat_d += 1;
                }
            }
            if (feat_d > d) {
                d = feat_d;
            }
        }
        dist[i] = d;
    }
    auto q = std::min_element(dist.begin(), dist.end());
    return q - dist.begin();
}

std::pair<std::vector<std::vector<short>>, int> BinCombAlg(std::vector<boost::dynamic_bitset<>> &X0, std::vector<boost::dynamic_bitset<>> &X1, int level)
{
    if (level > 0 && level <= 8) {
        std::cout << level << std::endl;
    }
    std::vector<std::vector<short>> opt_cl;
    int size_pred = 0;
    if ((level < 0 && (level + 100 > r)) || level > r) {
        return std::pair<std::vector<std::vector<short>>, int> (opt_cl, size_pred);
    }
    int q = find_object(X0, X1);
    std::vector<std::pair<std::vector<boost::dynamic_bitset<>>, int>> objects;
    int n = X0[0].size();
    int m0 = X0.size();
    int m1 = X1.size();
    for (int j = 0; j < n; j++) {
        std::vector<boost::dynamic_bitset<>> tmp;
        for (int i = 0; i < m0; i++) {
            if (X0[i][j] != X1[q][j]) {
                tmp.push_back(X0[i]);
            }
        }
        if (tmp.size() != 0) {
            objects.push_back(std::pair<std::vector<boost::dynamic_bitset<>>, int> (tmp, j));
        }
    }
    sort(objects.begin(), objects.end(), comp);
    for (auto ob : objects) {
        if (opt_cl.size() > 500000) {
            break;
        }
        std::vector<std::vector<short>> cur_cl;
        int cur_size = ob.first.size();
        if (cur_size < size_pred) {
            break;
        }
        int feature_num = ob.second;
        int feature_value = ob.first[0][feature_num];
        std::vector<short> el_cur_cl(n, -1);
        el_cur_cl[feature_num] = feature_value;
        cur_cl.push_back(el_cur_cl);
        std::vector<boost::dynamic_bitset<>> new_X_1;
        for (int i = 0; i < m1; i++) {
            if (X1[i][feature_num] == feature_value) {
                new_X_1.push_back(X1[i]);
            }
        }
        if (new_X_1.size() != 0) {
            std::pair<std::vector<std::vector<short>>, int> result = BinCombAlg(ob.first, new_X_1, level + 1);
            cur_size = result.second;
            cur_cl = result.first;
            for (unsigned i = 0; i < result.first.size(); i++) {
                for (int j = 0; j < n; j++) {
                    cur_cl[i][j] += (1 + el_cur_cl[j]);
                }
            }
        }
        if (cur_size > size_pred) {
            opt_cl = cur_cl;
            size_pred = cur_size;
        } else if (cur_size == size_pred) {
            for (unsigned i = 0; i < cur_cl.size(); i++) {
                opt_cl.push_back(cur_cl[i]);
            }
        }
    }
    return std::pair<std::vector<std::vector<short>>, int> (opt_cl, size_pred);
}

std::vector<std::string> predict(std::vector<std::vector<short>> &cl0, std::vector<std::vector<short>> &cl1, std::vector<double> &w0, std::vector<double> &w1, std::vector<boost::dynamic_bitset<>> &objects, std::pair<std::string, std::string> classes, std::string mode)
{
    std::vector<std::string> results;
    for (unsigned i = 0; i < objects.size(); i++) {
        double k0 = 0;
        double k1 = 0;
        for (unsigned j = 0; j < cl0.size(); j++) {
            int diff = 0;
            for (unsigned l = 0; l < objects[0].size(); l++) {
                if (cl0[j][l] != -1 && cl0[j][l] != objects[i][l]) {
                    diff++;
                    if (diff > 2) {
                        break;
                    }
                }
            }
            if (mode == "strict" && diff == 0) {
                k0 += w0[j];
            }
            if (mode == "soft" && diff <= 2) {
                k0 += w0[j];
            }
        }
        for (unsigned j = 0; j < cl1.size(); j++) {
            int diff = 0;
            for (unsigned l = 0; l < objects[0].size(); l++) {
                if (cl1[j][l] != -1 && cl1[j][l] != objects[i][l]) {
                    diff++;
                    if (diff > 2) {
                        break;
                    }
                }
            }
            if (mode == "strict" && diff == 0) {
                k1 += w1[j];
            }
            if (mode == "soft" && diff <= 2) {
                k1 += w1[j];
            }
        }
        if (k0 > k1) {
            results.push_back(classes.first);
        } else if (k0 < k1) {
            results.push_back(classes.second);
        } else {
            results.push_back("-");
        }
    }
    return results;
}

std::vector<std::pair<double, double>> predict_proba(std::vector<std::vector<short>> &cl0, std::vector<std::vector<short>> &cl1, std::vector<double> &w0, std::vector<double> &w1, std::vector<boost::dynamic_bitset<>> &objects)
{
    std::vector<std::pair<double, double>> results;
    for (unsigned i = 0; i < objects.size(); i++) {
        double k0 = 0;
        double k1 = 0;
        for (unsigned j = 0; j < cl0.size(); j++) {
            int diff = 0;
            for (unsigned l = 0; l < objects[0].size(); l++) {
                if (cl0[j][l] != -1 && cl0[j][l] != objects[i][l]) {
                    diff++;
                    if (diff > 2) {
                        break;
                    }
                }
            }
            if (diff == 0) {
                k0 += w0[j];
            }
        }
        for (unsigned j = 0; j < cl1.size(); j++) {
            int diff = 0;
            for (unsigned l = 0; l < objects[0].size(); l++) {
                if (cl1[j][l] != -1 && cl1[j][l] != objects[i][l]) {
                    diff++;
                    if (diff > 2) {
                        break;
                    }
                }
            }
            if (diff == 0) {
                k1 += w1[j];
            }
        }
        double s = k0 + k1;
        if (s < 1) {
            s = 2;
            k0 = 1;
            k1 = 1;
        }
        results.push_back({k0 / s, k1 / s});
    }
    
    return results;
}

double accuracy(std::vector<std::string> pred, std::vector<std::string> y_true) {
    double k = 0;
    for (unsigned i = 0; i < pred.size(); i++) {
        if (pred[i] == y_true[i]) {
            k++;
        }
    }
    return k / pred.size();
}

double accuracy1(std::vector<std::string>pred, std::vector<std::string> y_true)
{
    double k = 0;
    int n = 0;
    for (unsigned i = 0; i < pred.size(); i++) {
        if (pred[i] != "-") {
            if (pred[i] == y_true[i]) {
                k++;
            }
            n++;
        }
    }
    return k / n;
}

int main(int argc, char **argv){
    if (argc > 4) {
        r = 0;
        for (auto i: std::string(argv[4])) {
            r = r * 10 + (i - '0');
        }
    }
    if (argc > 1 && (std::string(argv[1]) == "rows" || std::string(argv[1]) == "cols")) {
        std::ofstream out("work_time_" + std::string(argv[1]) + "_simple_bit");
        for (int i = 0; i < 70; i++) {
            double work_time = 0;
            int p = 10;
            for (int k = 0; k < 10; k++) {
                std::vector<boost::dynamic_bitset<>> X0;
                std::vector<boost::dynamic_bitset<>> X1;
                std::string name_in = std::string(argv[1]) + "/" + std::string(argv[1]) + std::to_string(i + 1) + std::to_string(k);
                std::ifstream in(name_in);
                std::string line;
                if (in.is_open()) {
                    while (getline(in, line)) {
                        boost::dynamic_bitset<> str;
                        for (unsigned j = 0; j < line.size() - 1; j++) {
                            if (line[j] == '0') {
                                str.push_back(0);
                            } else if (line[j] == '1') {
                                str.push_back(1);
                            }
                        }
                        if (line[line.size() - 1] == '0') {
                            X0.push_back(str);
                        } else {
                            X1.push_back(str);
                        }
                    }
                    in.close();
                }
                if (X0.size() == 0 || X1.size() == 0) {
                    p -= 1;
                    continue;
                }
                double start_time = clock();
                if (argc > 2 && std::string(argv[2]) == "show") {
                    std::pair<std::vector<std::vector<short>>, int> result = BinCombAlg(X0, X1, 1);
                    std::cout << "|cov0| = " << result.second << std::endl;
                    result =  BinCombAlg(X1, X0, 1);
                    std::cout << "|cov1| = " << result.second << std::endl;
                } else {
                    std::pair<std::vector<std::vector<short>>, int> result =  BinCombAlg(X0, X1, -100);
                    std::cout << "|cov0| = " << result.second << std::endl;
                    result =  BinCombAlg(X1, X0, -100);
                    std::cout << "|cov1| = " << result.second << std::endl;
                }
                work_time += (clock() - start_time) / CLOCKS_PER_SEC;
            }
            std::cout << p << std::endl;
            std::cout << work_time << std::endl;
            out << work_time / p << std::endl;
        }
        out.close();
    } else if (argc > 1) {
        std::ofstream out(std::string(argv[1]) + "_time_simple_bit");
        std::vector<boost::dynamic_bitset<>> X0;
        std::vector<boost::dynamic_bitset<>> X1;
        std::pair<std::string, std::string> y("-", "-");
        std::ifstream in(argv[1]);
        std::string line;
        if (in.is_open()) {
            while (getline(in, line)) {
                boost::dynamic_bitset<> str;
                unsigned y_start = line.size() - 1;
                for (unsigned j = 0; j < line.size() - 1; j++) {
                    if (line[j] == '0') {
                        str.push_back(0);
                    } else if (line[j] == '1') {
                        str.push_back(1);
                    } else if (line[j] != ' ') {
                        y_start = j;
                        break;
                    }
                }
                std::string cur_y = line.substr(y_start, line.size() - y_start);
                if (y.first == "-") {
                    y.first = cur_y;
                    X0.push_back(str);
                } else if (y.first != cur_y && y.second == "-") {
                    y.second = cur_y;
                    X1.push_back(str);
                } else if (cur_y == y.first) {
                    X0.push_back(str);
                } else {
                    X1.push_back(str);
                }
            }
            in.close();
        }
        time_t start;
        struct tm *timeinfo_start;
        time(&start);
        timeinfo_start = localtime(&start);
        std::cout << "START: " << asctime(timeinfo_start) << std::endl;
        double start_time =  clock();
        std::vector<std::vector<short>> cl0;
        std::vector<std::vector<short>> cl1;
        std::vector<boost::dynamic_bitset<>> old_X0 = X0;
        std::vector<boost::dynamic_bitset<>> old_X1 = X1;
        std::vector<double> w0;
        std::vector<double> w1;
        std::cout << "cur sizes : " << X0.size() << '*' << X0[0].size() << ", " << X1.size() << '*' << X1[0].size() << std::endl;
        if (argc > 2 && std::string(argv[2]) == "show") {
            std::pair<std::vector<std::vector<short>>, int> result =  BinCombAlg(X0, X1, 1);
            std::set<std::vector<short>> unique(result.first.begin(), result.first.end());
            cl0.insert(cl0.end(), unique.begin(), unique.end());
            for (unsigned l = 0; l < unique.size(); l++) {
                w0.push_back(result.second);
            }
            std::cout << "|cov0| = " << result.second << std::endl;
            result =  BinCombAlg(X1, X0, 1);
            std::cout << "|cov1| = " << result.second << std::endl;
            std::set<std::vector<short>> unique1(result.first.begin(), result.first.end());
            cl1.insert(cl1.end(), unique1.begin(), unique1.end());
            for (unsigned l = 0; l < unique1.size(); l++) {
                w1.push_back(result.second);
            }
        } else {
            std::pair<std::vector<std::vector<short>>, int> result =  BinCombAlg(X0, X1, -100);
            std::cout << "|cov0| = " << result.second << std::endl;
            std::set<std::vector<short>> unique(result.first.begin(), result.first.end());
            cl0.insert(cl0.end(), unique.begin(), unique.end());
            for (unsigned l = 0; l < unique.size(); l++) {
                w0.push_back(result.second);
            }
            result =  BinCombAlg(X1, X0, -100);
            std::cout << "|cov1| = " << result.second << std::endl;
            std::set<std::vector<short>> unique1(result.first.begin(), result.first.end());
            cl1.insert(cl1.end(), unique1.begin(), unique1.end());
            for (unsigned l = 0; l < unique1.size(); l++) {
                w1.push_back(result.second);
            }
        }
        time_t end;
        struct tm *timeinfo_end;
        time(&end);
        timeinfo_end = localtime(&end);
        double work_time =  (clock() - start_time) / CLOCKS_PER_SEC;
        out << work_time << std::endl;
        out.close();
        std::cout << "FIT TIME: " << work_time << std::endl;
        std::cout << "END: " << asctime(timeinfo_end) << std::endl;
        std::vector<boost::dynamic_bitset<>> test;
        std::vector<std::string> y_true;
        std::ifstream in_test(argv[3]);
        if (in_test.is_open()) {
            while(getline(in_test, line)) {
                boost::dynamic_bitset<> str;
                unsigned y_start = line.size() - 1;
                for (unsigned j = 0; j < line.size() - 1; j++) {
                    if (line[j] == '0') {
                        str.push_back(0);
                    } else if (line[j] == '1') {
                        str.push_back(1);
                    } else if (line[j] != ' ') {
                        y_start = j;
                        break;
                    }
                }
                std::string cur_y = line.substr(y_start, line.size() - y_start);
                test.push_back(str);
                y_true.push_back(cur_y);
            }
            in_test.close();
        }
        start_time = clock();
        std::vector<std::string> y_pred = predict(cl0, cl1, w0, w1, test, y, "strict");
        std::cout << "before proba\n";
        std::vector<std::pair<double, double>> y_proba = predict_proba(cl0, cl1, w0, w1, test);
        work_time = (clock() - start_time) / CLOCKS_PER_SEC;
        std::cout << "PREDICT TIME: " << work_time << std::endl;
        double acc = accuracy(y_pred, y_true);
        double acc1 = accuracy1(y_pred, y_true);
        std::cout << "ACC: " << acc << std::endl;
        std::cout << "ACC1: " << acc1 << std::endl;
        std::ofstream out_pred(std::string(argv[1]) + "_pred_simple_bit");
        for (unsigned i = 0; i < y_pred.size(); i++) {
            out_pred << y_true[i] << " " << y_pred[i] << std::endl;
        }
        out_pred.close();

        std::ofstream out_proba(std::string(argv[1]) + "_proba_simple_bit");
        out_proba << y.first << " " << y.second << std::endl;
        for (unsigned i = 0; i < y_proba.size(); i++) {
            out_proba << y_true[i] << " " <<y_proba[i].first << " " << y_proba[i].second << std::endl;
        }
        out_proba.close();
    }
    return 0;
}
