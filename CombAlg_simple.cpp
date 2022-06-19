#include <iostream>
#include <vector>
#include <utility>
#include <algorithm>
#include <fstream>
#include <ctime>

bool comp(std::pair<std::vector<std::vector<short>>, int> &pat1, std::pair<std::vector<std::vector<short>>, int> &pat2)
{
    return pat1.first.size() > pat2.first.size();
}

int find_object(std::vector<std::vector<short>> &X0, std::vector<std::vector<short>> &X1)
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

std::pair<std::vector<std::vector<short>>, int> BinCombAlg(std::vector<std::vector<short>> &X0, std::vector<std::vector<short>> &X1, int level)
{
    if (level > 0 && level <= 8) {
        std::cout << level << std::endl;
    }
    std::vector<std::vector<short>> opt_cl;
    int size_pred = 0;
    int q = find_object(X0, X1);
    std::vector<std::pair<std::vector<std::vector<short>>, int>> objects;
    int n = X0[0].size();
    int m0 = X0.size();
    int m1 = X1.size();
    for (int j = 0; j < n; j++) {
        std::vector<std::vector<short>> tmp;
        for (int i = 0; i < m0; i++) {
            if (X0[i][j] != X1[q][j]) {
                tmp.push_back(X0[i]);
            }
        }
        if (tmp.size() != 0) {
            objects.push_back(std::pair<std::vector<std::vector<short>>, int> (tmp, j));
        }
    }
    sort(objects.begin(), objects.end(), comp);
    for (auto ob : objects) {
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
        std::vector<std::vector<short>> new_X_1;
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

std::vector<std::string> predict(std::vector<std::vector<short>> &cl0, std::vector<std::vector<short>> &cl1, std::vector<double> &w0, std::vector<double> &w1, std::vector<std::vector<short>> &objects, std::pair<std::string, std::string> classes, std::string mode)
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
        //k0 /= cl0.size();
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
    if (argc > 1 && (std::string(argv[1]) == "rows" || std::string(argv[1]) == "cols")) {
        std::ofstream out("work_time_" + std::string(argv[1]));
        for (int i = 0; i < 40; i++) {
            double work_time = 0;
            int p = 10;
            for (int k = 0; k < 10; k++) {
                std::vector<std::vector<short>> X0;
                std::vector<std::vector<short>> X1;
                std::string name_in = std::string(argv[1]) + std::to_string(i + 1) + std::to_string(k);
                std::ifstream in(name_in);
                std::string line;
                if (in.is_open()) {
                    while (getline(in, line)) {
                        std::vector<short> str;
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
                    BinCombAlg(X0, X1, 1);
                } else {
                    BinCombAlg(X0, X1, -100);
                }
                work_time += (clock() - start_time) / CLOCKS_PER_SEC;
            }
	    std::cout << p << std::endl;
	    std::cout << work_time << std::endl;
            out << work_time / p << std::endl;
        }
        out.close();
    } else if (argc > 1) {
        std::ofstream out(std::string(argv[1]) + "_time_simple");
        std::vector<std::vector<short>> X0;
        std::vector<std::vector<short>> X1;
	std::pair<std::string, std::string> y("-", "-");
        std::ifstream in(argv[1]);
        std::string line;
        if (in.is_open()) {
            while (getline(in, line)) {
		std::vector<short> str;
		unsigned y_start = 0;
                for (unsigned j = 0; j < line.size(); j++) {
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
	std::vector<std::vector<short>> old_X0 = X0;
	std::vector<std::vector<short>> old_X1 = X1;
	std::vector<double> w0;
	std::vector<double> w1;
	std::ofstream cl_x0(std::string(argv[1]) + "_" + y.first + "simple");
	std::ofstream cl_x1(std::string(argv[1]) + "_" + y.second + "simple");
        if (argc > 2 && std::string(argv[2]) == "show") {
            std::pair<std::vector<std::vector<short>>, int> result =  BinCombAlg(X0, X1, 1);
	    cl0.insert(cl0.end(), result.first.begin(), result.first.end());
	    cl_x0 << result.second << '\n';
	    for (unsigned l = 0; l < result.first.size(); l++) {
	        for (unsigned ll = 0; ll < result.first[0].size(); ll++) {
		    cl_x0 << result.first[l][ll] << " ";
	        }
		cl_x0 << '\n';
		w0.push_back(result.second);
	    }
	    result =  BinCombAlg(X1, X0, 1);
	    cl1.insert(cl1.end(), result.first.begin(), result.first.end());
	    cl_x1 << result.second << '\n';
            for (unsigned l = 0; l < result.first.size(); l++) {
                for (unsigned ll = 0; ll < result.first[0].size(); ll++) {
                    cl_x1 << result.first[l][ll] << ' ';
                }
                cl_x1 << '\n';
                w1.push_back(result.second);
            }
        } else {
            std::pair<std::vector<std::vector<short>>, int> result =  BinCombAlg(X0, X1, -100);
            cl0.insert(cl0.end(), result.first.begin(), result.first.end());
            cl_x0 << result.second << '\n';
            for (unsigned l = 0; l < result.first.size(); l++) {
                for (unsigned ll = 0; ll < result.first[0].size(); ll++) {
                    cl_x0 << result.first[l][ll] << ' ';
                }
                cl_x0 << '\n';
                w0.push_back(result.second);
            }
            result =  BinCombAlg(X1, X0, -100);
            cl1.insert(cl1.end(), result.first.begin(), result.first.end());
            cl_x1 << result.second << '\n';
            for (unsigned l = 0; l < result.first.size(); l++) {
                for (unsigned ll = 0; ll < result.first[0].size(); ll++) {
                    cl_x1 << result.first[l][ll] << ' ';
                }
                cl_x1 << '\n';
                w1.push_back(result.second);
            }
        }
	std::cout << "FOUND\n";
	cl_x0.close();
	cl_x1.close();
	time_t end;
        struct tm *timeinfo_end;
        time(&end);
        timeinfo_end = localtime(&end);
        double work_time =  (clock() - start_time) / CLOCKS_PER_SEC;
        out << work_time << std::endl;
        out.close();
	std::cout << "FIT TIME: " << work_time << std::endl;
	std::cout << "END: " << asctime(timeinfo_end) << std::endl;
	std::vector<std::vector<short>> test;
        std::vector<std::string> y_true;
	std::ifstream in_test(argv[3]);
	if (in_test.is_open()) {
            while(getline(in_test, line)) {
                std::vector<short> str;
		unsigned y_start = 0;
                for (unsigned j = 0; j < line.size(); j++) {
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
        work_time = (clock() - start_time) / CLOCKS_PER_SEC;
	std::cout << "PREDICT TIME: " << work_time << std::endl;
        std::ofstream pred_str(std::string(argv[3]) + "_predict_simple");
        if (pred_str.is_open()) {
            for (unsigned i = 0; i < y_pred.size(); i++) {
                pred_str << y_pred[i] << std::endl;
            }
        }
        pred_str.close();
        double acc = accuracy(y_pred, y_true);
        double acc1 = accuracy1(y_pred, y_true);
        std::cout << "ACC: " << acc << std::endl;
        std::cout << "ACC1: " << acc1 << std::endl;
        std::ofstream out_acc(std::string(argv[1]) + "_predict_acc_simple");
        out_acc << acc << '\n' << acc1 << std::endl;
        out_acc.close();
    }
    return 0;
}
