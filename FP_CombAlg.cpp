#include <iostream>
#include <vector>
#include <utility>
#include <algorithm>
#include <fstream>
#include <ctime>
#include <map>
#include <boost/dynamic_bitset.hpp>


unsigned all_new = 0;
unsigned all_delete = 0;
bool glob_flag = false;
bool x1_flag = false;
unsigned x0_size = 0;
int max_level = 0;

class Node
{
public:
    Node *prev;
    std::vector<Node *> next;
    int feature_num;
    unsigned val;
    Node(): prev(nullptr), next(0), feature_num(-1), val(0) {}
    Node(const Node&v): prev(v.prev), next(v.next), feature_num(v.feature_num), val(v.val) {}
    const Node & operator=(const Node& v) {
        prev = v.prev;
        next = v.next;
        feature_num = v.feature_num;
        val = v.val;
        return *this;
    }
};

class MainNode
{
public:
    int feature_num;
    unsigned counter;
    std::vector<Node *> node;
    MainNode(): feature_num(-1), counter(0), node(0) {}
};

bool comp1(MainNode &n1, MainNode &n2)
{
    return n1.counter > n2.counter;
}

bool comp2(std::pair<int, int> &n1, std::pair<int, int> &n2)
{
    if (n1.first < n2.first) {
        return true;
    }
    if (n1.first == n2.first) {
        return n1.second > n2.second;
    }
    return false;
}


void show_root(Node *root)
{
    std::cout << "+++++++++++++++++++++\n";
    for (auto ob: root->next) {
        std::cout << ob << " " << ob->feature_num << std::endl;
    }
    std::cout << "+++++++++++++++++++++\n";
}

void show_X(std::vector<boost::dynamic_bitset<>> &X, std::vector<unsigned> cur_list)
{
    std::cout << "===================\n";
    for (auto i: cur_list) {
        for (unsigned j = 1; j < X[i].size(); j += 2) {
            std::cout << X[i][j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << "===================\n";
}

void show_bound(unsigned bound, std::vector<boost::dynamic_bitset<>> &X1)
{
    for (unsigned i = 1; i < X1[0].size(); i += 2) {
        std::cout << X1[bound][i] << " ";
    }
    std::cout << std::endl;
}

void show_tree(Node *root)
{
    for (auto ob: root->next) {
        std::cout << "adr: " << ob << ", prev: " << root << ", feat_num: " << ob->feature_num << ", count: " << ob->val << std::endl;
        show_tree(ob);
    }
}

void show_mainnodes(std::vector<MainNode> &m)
{
    for (auto v: m) {
        std::cout << "feature_num: " << v.feature_num << ", counter = " << v.counter << std::endl;
    }
}


int find_object(std::vector<boost::dynamic_bitset<>> &X0, std::vector<boost::dynamic_bitset<>> &X1, std::vector<unsigned> &cur_list1, std::vector<MainNode> &mainnodes) //, std::map<int, int> &nodes_map)
{
    int n = mainnodes.size();
    std::vector<std::pair<int, int>> dist;
    for (auto i: cur_list1) {
        int d = 0;
        int nodes_num = 0;
        for (int j = 0; j < n; j++) {
            int feat_d = 0;
            int cur_num = mainnodes[j].feature_num;
            if (X1[i][cur_num] == 0) {
                feat_d = mainnodes[j].counter;
            }
            if (feat_d >= d) {
                d = feat_d;
                nodes_num = j;
            }
        }
        dist.push_back(std::pair<int, int> {d, nodes_num});
    }
    auto q = std::min_element(dist.begin(), dist.end(), comp2);
    return q - dist.begin();
}

Node *build(std::vector<boost::dynamic_bitset<>> &X0, std::vector<MainNode> &mainnodes)
{
    Node *root = new Node; // главный корень, не соответствующий никакому признаку
    all_new++; // выделение памяти
    for (auto ob : X0) { // проходим по всем прецедентам класса К
        auto cur_node = root;
        for (unsigned i = 0; i < mainnodes.size(); i++) { // сохраняем информацию об объекте
            if (ob[mainnodes[i].feature_num] == 0) { // в данном объекте нет такого признака
                continue;
            }
            bool flag = true;
            for (unsigned j = 0; j < cur_node->next.size(); j++) { // вершина уже есть
                if (cur_node->next[j]->feature_num == mainnodes[i].feature_num) { 
                    flag = false; // новую вершину создавать не нужно
                    cur_node = cur_node->next[j];
                    cur_node->val++; // обновляем значение счетчика
                    break;
                }
            }
            if (flag) { // вершины нет, создаем новую
                Node *nn = new Node;
                all_new++; // выделение памяти
                nn->prev = cur_node; // обозначаем родителя
                nn->feature_num = mainnodes[i].feature_num; // номер признака
                nn->val = 1; // счетчик
                mainnodes[i].node.push_back(nn); // записываем в список вершин соответствующего признака
                cur_node->next.push_back(nn); // записываем в список вершин
                cur_node = nn; // спускаемся вниз
            }
       }
   }
   return root; // возвращаем корень
}

std::pair<Node *, std::vector<MainNode>> cond_tree(std::vector<MainNode> &mainnodes, Node *old_root, unsigned i)
{
    std::vector<MainNode> new_mainnodes(i); // создаем новый список признаков
    if (new_mainnodes.size() == 0) { // исключаем случай корня
        return std::pair<Node *, std::vector<MainNode>>(nullptr, new_mainnodes);
    }
    std::copy(mainnodes.begin(), mainnodes.begin() + i, new_mainnodes.begin()); // копируем старый список признаков в новый
    for (unsigned j = 0; j < new_mainnodes.size(); j++) {
        new_mainnodes[j].node.clear(); // стираем ссылки, потому что они будут новые
        new_mainnodes[j].counter = 0; // обнуляем счетчик
    }
    std::map<int, unsigned> nums_indexes; // отображение реальных номеров признаков в текущие индексы
    for (unsigned j = 0; j < new_mainnodes.size(); j++) {
        nums_indexes.insert({new_mainnodes[j].feature_num, j});
    }
    std::map<Node *, Node *> old_new_adresses; // отображение старых адресов в новые, чтобы не раздваивать вершины
    Node *root = new Node; // создаем новый корень
    old_new_adresses[old_root] = root;
    all_new++; // выделение памяти
    for (auto p : mainnodes[i].node) { // обход снизу вверх
        if (old_new_adresses.find(p->prev) != old_new_adresses.end()) { // такую вершину мы уже создали, обновляем счетчик
            Node *tmp = old_new_adresses[p->prev]; // адрес новой вершины сверху
            while (tmp->prev != nullptr and tmp->prev->feature_num != -1) { // пока не уткнулись в корень
                tmp->val += p->val;
                tmp = tmp->prev;
            }
            tmp->val += p->val;
            continue;
        }
        if (p->prev->feature_num == -1) { // если уткнулись в корень
            continue;
        }
        Node *cur_new = new Node;
        all_new++; // выделение памяти
        old_new_adresses[p->prev] = cur_new; // отображение старого родителя в нового родителя
        Node *cur_old = p->prev; // поднялись в родителя
        cur_new->val = p->val; // заполняем счетчик
        cur_new->feature_num = cur_old->feature_num; // заполняем номер признака
        unsigned ind = nums_indexes[cur_new->feature_num]; // получаем индекс мэйнноды
        new_mainnodes[ind].node.push_back(cur_new); // пихаем в нее созданную вершину
        while (cur_old->prev->feature_num != -1) { // пока не уткнулись в корень
            if (old_new_adresses.find(cur_old->prev) != old_new_adresses.end()) { // если такая вершина уже есть
                cur_new->prev = old_new_adresses[cur_old->prev]; // новый родитель
                cur_old = cur_old->prev; // старый родитель
                if (std::find(cur_new->prev->next.begin(), cur_new->prev->next.end(), cur_new) == cur_new->prev->next.end()) { // если в новом родителе нет текущего ребенка
                    cur_new->prev->next.push_back(cur_new);
                }
                Node *tmp = cur_new->prev;
                while (tmp->prev != nullptr and tmp->prev->feature_num != -1) { // обновляем счетчики на текущем пути
                    tmp->val += p->val;
                    tmp = tmp->prev;
                }
                tmp->val += p->val;
                cur_new = tmp;//cur_new->prev; // идем наверх
                break;
            } else {
                Node *new_node = new Node;
                all_new++; // выделение памяти
                old_new_adresses[cur_old->prev] = new_node; // заполняем адрес
                cur_new->prev = new_node; // заполняем все поля
                new_node->feature_num = cur_old->prev->feature_num;
                new_node->next.push_back(cur_new);
                cur_old = cur_old->prev;
                new_node->val = cur_new->val;
                cur_new = new_node;
                ind = nums_indexes[cur_new->feature_num];
                new_mainnodes[ind].node.push_back(cur_new);
            }
        }
        if (std::find(root->next.begin(), root->next.end(), cur_new) == root->next.end()) {
            if (cur_new->feature_num != -1) {
                root->next.push_back(cur_new);
                cur_new->prev = root;
            }
        }
    }
    auto new_end = std::remove_if(new_mainnodes.begin(), new_mainnodes.end(), [](MainNode m) { return m.node.size() == 0; });
    new_mainnodes.erase(new_end, new_mainnodes.end());
    for (unsigned j = 0; j < new_mainnodes.size(); j++) {
        for (auto p : new_mainnodes[j].node) {
            new_mainnodes[j].counter += p->val;
        }
    }
    return std::pair<Node *, std::vector<MainNode>>(root, new_mainnodes);
}

void del(Node *root)
{
    if (root == nullptr) {
        return;
    }
    if (root->next.size() == 0) {
        all_delete++;
        delete root;
    } else {
        for (auto p : root->next) {
            del(p);
        }
        all_delete++;
        delete root;
    }
}

std::pair<unsigned, std::vector<std::vector<int>>> find(std::vector<boost::dynamic_bitset<>> &X0, std::vector<boost::dynamic_bitset<>> &X1, Node *root, 
std::vector<MainNode> &mainnodes, std::vector<unsigned> cur_list1, unsigned maximum, int level, int prev_f, std::map<int, int> &nodes_map)
{
    unsigned cur_max = 0;
    std::vector<std::vector<int>> features;
    unsigned bound = cur_list1[find_object(X0, X1, cur_list1, mainnodes)];
    for (int i = 0; i < int(mainnodes.size()); i++) {
        if (X1[bound][mainnodes[i].feature_num] == 1) {
            continue;
        }
        if (mainnodes[i].counter >= 0.1 * X0.size() && mainnodes[i].counter >= maximum) {
            std::vector<unsigned> new_list1; // номера объектов из X1, которые пока не классифицировали
            for (auto x_num: cur_list1) {
                if (X1[x_num][mainnodes[i].feature_num] == 1) {
                    new_list1.push_back(x_num);
                }
            }
            if (new_list1.size() == 0) { // если таких объектов нет
                if (mainnodes[i].counter > maximum) {
                    features.clear();
                    maximum = mainnodes[i].counter;
                }
                cur_max = maximum;
                features.push_back(std::vector<int>{mainnodes[i].feature_num});
                continue; // переход к следующему признаку
            }
            std::pair<Node *, std::vector<MainNode>> new_pair = cond_tree(mainnodes, root, i); // строим условное дерево
            Node *new_tree = new_pair.first;
            std::vector<MainNode> new_mainnodes = new_pair.second;
            if (new_mainnodes.size() == 0) { // если это дерево пустое
                if (new_tree) {
                    del(new_tree);
                }
                continue;
            }
            auto res = find(X0, X1, new_tree, new_mainnodes, new_list1, maximum, level + 1, mainnodes[i].feature_num, nodes_map); // запускаем рекурсию
            unsigned count = res.first;
            std::vector<std::vector<int>> cur_features = res.second;
            del(new_tree);
            if (count == 0 || count < maximum) {
                continue;
            }
            if (count > maximum) {
                maximum = count;
                features.clear();
            }
            if (count >= maximum) {
                cur_max = maximum;
                for (auto feat : cur_features) {
                    feat.push_back(mainnodes[i].feature_num);
                    features.push_back(feat);
                }
            }
            if (features.size() > 500000) {
                return std::pair<unsigned, std::vector<std::vector<int>>>(cur_max, features);
            }
        }
    }
    return std::pair<unsigned, std::vector<std::vector<int>>>(cur_max, features);
}


std::pair<std::vector<std::vector<int>>, int> FP(std::vector<boost::dynamic_bitset<>> &X0, std::vector<boost::dynamic_bitset<>> &X1)
{
    std::vector<MainNode> mainnodes(X0[0].size());
    for (auto ob : X0) { // строим левые вершины, подсчитываем число признаков
        for (unsigned i = 0; i < X0[0].size(); i++) {
            mainnodes[i].counter += ob[i];
        }
    }
    for (unsigned i = 0; i < X0[0].size(); i++) {
        mainnodes[i].feature_num = i; // заполняем номер признака
    }
    x0_size = X0.size();
    auto new_end = std::remove_if(mainnodes.begin(), mainnodes.end(), [](MainNode m) { return m.counter < 0.1 * x0_size; });
    mainnodes.erase(new_end, mainnodes.end());
    sort(mainnodes.begin(), mainnodes.end(), comp1); // сортируем список левых вершин по убыванию
    std::map<int, int> nodes_map;
    for (unsigned i = 0; i < mainnodes.size(); i++) {
        nodes_map[mainnodes[i].feature_num] = i;
    }
    Node *root = build(X0, mainnodes); // строим основное дерево (диамически)
    unsigned maximum = 0;
    std::vector<unsigned> cur_list1;
    for (unsigned i = 0; i < X1.size(); i++) {
        cur_list1.push_back(i);
    }
    auto res = find(X0, X1, root, mainnodes, cur_list1, maximum, 0, -1, nodes_map);
    std::vector<std::vector<int>> features = res.second;
    unsigned count = res.first;
    int c = count;
    del(root);
    return std::pair<std::vector<std::vector<int>>, int>(features, c);
}

std::vector<std::string> predictFP(std::vector<std::vector<int>> &cl0, std::vector<std::vector<int>> &cl1, std::vector<double> &w0, std::vector<double> &w1, std::vector<boost::dynamic_bitset<>> &objects, std::pair<std::string, std::string> classes)
{
    std::vector<std::string> results;
    for (unsigned i = 0; i < objects.size(); i++) {
        double k0 = 0;
        double k1 = 0;
        for (unsigned j = 0; j < cl0.size(); j++) {
            int diff = 0;
            for (auto feat_num: cl0[j]) {
                if (objects[i][feat_num] == 0) {
                    diff++;
                    break;
                }
            }
            if (diff == 0) {
                k0 += w0[j];
            }
        }
        for (unsigned j = 0; j < cl1.size(); j++) {
            int diff = 0;
            for (auto feat_num: cl1[j]) {
                if (objects[i][feat_num] == 0) {
                    diff++;
                    break;
                }
            }
            if (diff == 0) {
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

std::vector<std::pair<double, double>> predict_proba(std::vector<std::vector<int>> &cl0, std::vector<std::vector<int>> &cl1, std::vector<double> &w0, std::vector<double> &w1, std::vector<boost::dynamic_bitset<>> &objects)
{
    std::vector<std::pair<double, double>> results;
    for (unsigned i = 0; i < objects.size(); i++) {
        double k0 = 0;
        double k1 = 0;
        for (unsigned j = 0; j < cl0.size(); j++) {
            int diff = 0;
            for (auto feat_num: cl0[j]) {
                if (objects[i][feat_num] == 0) {
                    diff++;
                    break;
                }
            }
            if (diff == 0) {
                k0 += w0[j];
            }
        }
        for (unsigned j = 0; j < cl1.size(); j++) {
            int diff = 0;
            for (auto feat_num: cl1[j]) {
                if (objects[i][feat_num] == 0) {
                    diff++;
                    break;
                }
            }
            if (diff == 0) {
                k1 += w1[j];
            }
        }
        double s = k0 + k1;
        if (s < 0.5) {
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
    if (argc > 1 && (std::string(argv[1]) == "rows" || std::string(argv[1]) == "cols")) {
        std::ofstream out("work_time_new_" + std::string(argv[1]) + "_FP_bit_min");
        std::pair<std::string, std::string> y("0", "1");
        for (int i = 0; i < 70; i++) {
            std::cout << "##############################################\n";
            double work_time = 0;
            int p = 10;
            for (int k = 0; k < 10; k++) {
                std::cout << "*****************************************\n";
                all_new = 0;
                all_delete = 0;
                std::vector<boost::dynamic_bitset<>> X0;
                std::vector<boost::dynamic_bitset<>> X1;
                std::string name_in = "new_" + std::string(argv[1]) + '/' + std::string(argv[1]) + std::to_string(i + 1) + std::to_string(k);
                std::ifstream in(name_in);
                std::string line;
                if (in.is_open()) {
                    while (getline(in, line)) {
                        boost::dynamic_bitset<> str;
                        for (unsigned j = 0; j < line.size() - 1; j++) {
                            if (line[j] == '0') {
                                str.push_back(1);
                                str.push_back(0);
                            } else if (line[j] == '1') {
                                str.push_back(0);
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
                std::vector<std::vector<int>> cl0;
	            std::vector<std::vector<int>> cl1;
	            std::vector<boost::dynamic_bitset<>> old_X0 = X0;
	            std::vector<boost::dynamic_bitset<>> old_X1 = X1;
	            std::vector<double> w0;
	            std::vector<double> w1;
                double start_time = clock();
                int n_iter = 0;
                int cov0 = 10;
                int cov1 = 10;
                while (X0.size() > 0 && X1.size() > 0 && (X0.size() > 1 || X1.size() > 1) && k < 20 && (cov0 > 0 || cov1 > 0)) { 
                    std::cout << "cur sizes : " << X0.size() << '*' << X0[0].size() << ", " << X1.size() << '*' << X1[0].size() << std::endl;
                    n_iter++;
                    std::cout << "ITER:" << n_iter << std::endl;
                    max_level = 0;
                    std::pair<std::vector<std::vector<int>>, int> result;
                    if (X0.size() > 1 && cov0 > 0) {
                        result = FP(X0, old_X1);
                        if (all_new != all_delete) {
                            std::cout << "new/delete: " << all_new << "/" << all_delete << std::endl;
                            std::cout << "___________________________\n";
                            std::cout << "Memory Error!\n";
                            return 0;
                        }
                        all_new = 0;
                        all_delete = 0;
                        cl0.insert(cl0.end(), result.first.begin(), result.first.end());
                        cov0 = result.second;
                        std::cout << "|cov0| = " << result.second << std::endl;
                        for (unsigned l = 0; l < result.first.size(); l++) {
                            w0.push_back(result.second);
                        }
                    } else {
                        cov0 = 0;
                    }
                    if (X1.size() > 1 && cov1  > 0) {
                        result = FP(X1, old_X0);
                        cl1.insert(cl1.end(), result.first.begin(), result.first.end());
                        cov1 = result.second;
                        std::cout << "|cov1| = " << result.second << std::endl;
                        for (unsigned l = 0; l < result.first.size(); l++) {
                            w1.push_back(result.second);
                        }
                        
                    } else {
                        cov1 = 0;
                    }
                    std::vector<std::string> pred_X0 = predictFP(cl0, cl1, w0, w1, X0, y);
                    std::vector<boost::dynamic_bitset<>> new_X0;
                    for (unsigned l = 0; l < pred_X0.size(); l++) {
                        if (pred_X0[l] == "-") {
                            new_X0.push_back(X0[l]);
                        }
                    }
	                X0 = new_X0;
	                std::vector<std::string> pred_X1 = predictFP(cl0, cl1, w0, w1, X1, y);
                    if (all_new != all_delete) {
                        std::cout << "new/delete: " << all_new << "/" << all_delete << std::endl;
                        std::cout << "___________________________\n";
                        std::cout << "Memory Error!\n";
                        return 0;
                    }
                    all_new = 0;
                    all_delete = 0;
                    std::vector<boost::dynamic_bitset<>> new_X1;
                    for (unsigned l = 0; l < pred_X1.size(); l++) {
                        if (pred_X1[l] == "-") {
                            new_X1.push_back(X1[l]);
                        }
                    }
                    X1 = new_X1;
                }
                work_time += (clock() - start_time) / CLOCKS_PER_SEC;
            }
	        std::cout << p << std::endl;
	        std::cout << work_time << std::endl;
            out << work_time / p << std::endl;
        }
        out.close();
    } else if (argc > 1) {
        std::ofstream out(std::string(argv[1]) + "_time_FP_bit_min");
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
                        str.push_back(1);
                        str.push_back(0);
                    } else if (line[j] == '1') {
                        str.push_back(0);
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
        std::vector<std::vector<int>> cl0;
        std::vector<std::vector<int>> cl1;
        std::vector<boost::dynamic_bitset<>> old_X0 = X0;
        std::vector<boost::dynamic_bitset<>> old_X1 = X1;
        std::vector<unsigned> num_X0;
        std::vector<unsigned> num_X1;
        for (unsigned i = 0; i < X0.size(); i++) {
            num_X0.push_back(i);
        }
        for (unsigned i = 0; i < X1.size(); i++) {
            num_X1.push_back(i);
        }
        std::vector<double> w0;
        std::vector<double> w1;
        int k = 0;
        int cov0 = 10;
        int cov1 = 10;
        std::ofstream proc(std::string(argv[1]) + "_process_FP_bit_min");
        std::ofstream inf(std::string(argv[1]) + "_inf_FP_bit_min");
        inf << "number for K | rg for K | number for !K | rg for !K\n";
        while ((X0.size() > 0 || X1.size() > 0) && (X0.size() > 6 || X1.size() > 6) && k < 20 && (cov0 > 0 || cov1 > 0)) {
            
            std::cout << "cur sizes : " << X0.size() << '*' << X0[0].size() << ", " << X1.size() << '*' << X1[0].size() << std::endl;
            k++;
            std::pair<std::vector<std::vector<int>>, int> result;
            std::cout << "ITER:" << k << std::endl;
            if (X0.size() > 6 && cov0 > 0) {
                result =  FP(X0, old_X1);
                cov0 = result.second;
                if (all_new != all_delete) {
                    std::cout << "new/delete: " << all_new << "/" << all_delete << std::endl;
                    std::cout << "___________________________\n";
                    std::cout << "Memory Error!\n";
                    return 0;
                }
                all_new = 0;
                all_delete = 0;

                cl0.insert(cl0.end(), result.first.begin(), result.first.end());
                std::cout << "|cov0| = " << result.second << std::endl;
                cov0 = result.second;
                for (unsigned l = 0; l < result.first.size(); l++) {
                    w0.push_back(result.second);
                }
                max_level = 0;
                int number = result.first.size();
                double rg = 0;
                for (auto r: result.first) {
                    rg += r.size();
                }
                rg /= number;
                inf << number << " | " << rg << " | ";
            } else {
                cov0 = 0;
                inf << " - | - | ";
            }
            if (X1.size() > 6 && cov1 > 0) {
                result =  FP(X1, old_X0);
                cov1 = result.second;
                if (all_new != all_delete) {
                    std::cout << "new/delete: " << all_new << "/" << all_delete << std::endl;
                    std::cout << "___________________________\n";
                    std::cout << "Memory Error!\n";
                    return 0;
                }
                all_new = 0;
                all_delete = 0;

                cl1.insert(cl1.end(), result.first.begin(), result.first.end());
                std::cout << "|cov1| = " << result.second << std::endl;
                for (unsigned l = 0; l < result.first.size(); l++) {
                    w1.push_back(result.second);
                }
                cov1 = result.second;
                max_level = 0;
                int number1 = result.first.size();
                double rg1 = 0;
                for (auto r: result.first) {
                    rg1 += r.size();
                }
                rg1 /= number1;
                inf << number1 << " | " << rg1 << "\n";
            } else {
                cov1 = 0;
                inf << " - | -\n";;
            }
            proc << "i " << k << std::endl;
	        std::vector<std::string> pred_X0 = predictFP(cl0, cl1, w0, w1, X0, y);
	        std::vector<boost::dynamic_bitset<>> new_X0;
	        for (unsigned i = 0; i < pred_X0.size(); i++) {
	            if (pred_X0[i] == "-") {
		            new_X0.push_back(X0[i]);
		        } else {
                    proc << num_X0[i] << " ";
                }
	        }
            proc << std::endl;
            for (unsigned i = 0, t = 0; i < pred_X0.size(); i++) {
	            if (pred_X0[i] != "-") {
                    num_X0.erase(num_X0.begin() + (i - t));
                    t++;
                }
	        }
	        X0 = new_X0;

	        std::vector<std::string> pred_X1 = predictFP(cl0, cl1, w0, w1, X1, y);
            std::vector<boost::dynamic_bitset<>> new_X1;
            for (unsigned i = 0; i < pred_X1.size(); i++) {
                if (pred_X1[i] == "-") {
                    new_X1.push_back(X1[i]);
                } else {
                    proc << num_X1[i] << " ";
                }
	        }
            for (unsigned i = 0, t = 0; i < pred_X1.size(); i++) {
	            if (pred_X1[i] != "-") {
                    num_X1.erase(num_X1.begin() + (i - t));
                    t++;
                }
	        }
            proc << std::endl;
            X1 = new_X1;
        }
        std::cout << "not classified : " << X0.size() << '*' << X0[0].size() << ", " << X1.size() << '*' << X1[0].size() << std::endl;
        time_t end;
        struct tm *timeinfo_end;
        time(&end);
        timeinfo_end = localtime(&end);
        double work_time =  (clock() - start_time) / CLOCKS_PER_SEC;
        out << work_time << std::endl;
        out.close();
        proc.close();
        std::cout << "FIT TIME: " << work_time << std::endl;
        std::cout << "END: " << asctime(timeinfo_end) << std::endl;
        std::vector<boost::dynamic_bitset<>> test;
        std::vector<std::string> y_true;
        std::ifstream in_test(argv[3]);
        if (in_test.is_open()) {
            while(getline(in_test, line)) {
                boost::dynamic_bitset<> str;
                unsigned y_start = line.size() - 1;
                for (unsigned j = 0; j < line.size(); j++) {
                    if (line[j] == '0') {
                        str.push_back(1);
                        str.push_back(0);
                    } else if (line[j] == '1') {
                        str.push_back(0);
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
        std::vector<std::string> y_pred = predictFP(cl0, cl1, w0, w1, test, y);
        std::vector<std::pair<double, double>> y_proba = predict_proba(cl0, cl1, w0, w1, test);
        work_time = (clock() - start_time) / CLOCKS_PER_SEC;
        std::cout << "PREDICT TIME: " << work_time << std::endl;
        double acc = accuracy(y_pred, y_true);
        double acc1 = accuracy1(y_pred, y_true);
        std::cout << "ACC: " << acc << std::endl;
        std::cout << "ACC1: " << acc1 << std::endl;
        std::cout << "k iterations: " << k << std::endl;
        std::ofstream out_pred(std::string(argv[1]) + "_pred_FP_bit_min");
        for (unsigned i = 0; i < y_pred.size(); i++) {
            out_pred << y_true[i] << " " << y_pred[i] << std::endl;
        }
        out_pred.close();

        std::ofstream out_proba(std::string(argv[1]) + "_proba_FP_bit_min");
        out_proba << y.first << " " << y.second << std::endl;
        for (unsigned i = 0; i < y_proba.size(); i++) {
            out_proba << y_true[i] << " " <<y_proba[i].first << " " << y_proba[i].second << std::endl;
        }
        out_proba.close();
    }
    return 0;
}
