#include<unordered_map>
#include<stack>
#include<vector>
#include<string>
#include<iostream>
#include<memory>
#include<time.h>
#include<random>
#include<algorithm>
#include<cassert>

inline int rand_between(int a, int b) {
    if (a >= b) {
        return a;
    }

    return rand() % (b - a) + a;
}

// def rand_parts(seq, n, l):
//     indices = xrange(len(seq) - (l - 1) * n)
//     result = []
//     offset = 0
//     for i in sorted(random.sample(indices, n)):
//         i += offset
//         result.append(seq[i:i+l])
//         offset += l - 1
//     return result

std::string generate_random_string(int length) {
    const char BASES[4] = {'A', 'C', 'G', 'T'};
    std::string result = "";

    for (int i = 0; i < length; i++) {
        result += BASES[rand_between(0, 4)];
    }

    return result;
}

f(int a, int b) {
    return a < b;
}

std::string add_repeats(std::string seq, std::string repeat, int n) {
    int l = repeat.length();
    int range_size = seq.length() - (l - 1) * n;
    std::vector<int> range;
    std::vector<int> rand_sample;
    
    for(int i = 0; i < seq.length() - (l - 1) * n; i++) {
        range.push_back(i);
    }
    //std::experimental::sample(range.begin(), range.end(), std::back_inserter(rand_sample), n, std::mt19937{std::random_device{}()});

    std::random_shuffle(range.begin(), range.end());
    std::vector<int>::iterator it = range.begin();
    for (int i = 0; i < n; i++) {
        rand_sample.push_back(*it);
        it++;
    }
    std::sort(rand_sample.begin(), rand_sample.end(), f);

    std::string result = seq;
    int offset = 0;
    for (std::vector<int>::iterator it = rand_sample.begin(); it != rand_sample.end(); it++) {
        result = result.substr(0, *it + offset) + repeat + result.substr(*it + offset + l);
        offset += l - 1;
    }
    
    return result;
}

struct Node {
    Node(std::string value, std::unordered_map<std::shared_ptr<Node>, int> next, std::unordered_map<std::shared_ptr<Node>, int> prev) {
        this->value = value;
        this->next = next;
        this->prev = prev;
    }
    
    std::string value;
    std::unordered_map<std::shared_ptr<Node>, int> prev, next;
};

class Graph {
    public:
        Graph(std::string str, int read_length = 10) {
            std::string prefix, suffix;
            std::shared_ptr<Node> prefix_node, suffix_node;

            for (int i = 0; i < str.length() - read_length + 1; i++) {
                prefix = str.substr(i, read_length - 1);
                suffix = str.substr(i + 1, read_length - 1);
                prefix_node = add_vertice(prefix);
                suffix_node = add_vertice(suffix);
                add_edge(prefix_node, suffix_node);
            }
        }

        std::shared_ptr<Node> add_vertice(std::string value) {
            if (this->vertices.find(value) == this->vertices.end())
                this->vertices[value] = std::make_shared<Node>(value, std::unordered_map<std::shared_ptr<Node>, int>(), std::unordered_map<std::shared_ptr<Node>, int>());

            return this->vertices[value];
        }

        void add_edge(std::shared_ptr<Node> n1, std::shared_ptr<Node> n2) {
            if (n1->next.find(n2) == n1->next.end()) {
                n1->next[n2] = 1;
                n2->prev[n1] = 1;
            } else {
                n1->next[n2] += 1;
                n2->prev[n1] += 1;
            }
        }

        void remove_edge(std::shared_ptr<Node> n1, std::shared_ptr<Node> n2) {
            n1->next[n2] -= 1;
            n2->prev[n1] -= 1;

            if (n1->next[n2] <= 0) {
                n1->next.erase(n2);
                n2->prev.erase(n1);
            }
        }

        int calculate_indegree(std::shared_ptr<Node> node) {
            int sum = 0;
            for (std::unordered_map<std::shared_ptr<Node>, int>::iterator it = node->prev.begin(); it != node->prev.end(); it++) {
                sum += it->second;
            }

            return sum;
        }

        int calculate_outdegree(std::shared_ptr<Node> node) {
            int sum = 0;
            for (std::unordered_map<std::shared_ptr<Node>, int>::iterator it = node->next.begin(); it != node->next.end(); it++) {
                sum += it->second;
            }

            return sum;
        }


        std::vector<std::shared_ptr<Node>> find_eulerian_path() {
            std::stack<std::shared_ptr<Node>> node_stack;
            std::vector<std::shared_ptr<Node>> node_path;
            std::shared_ptr<Node> odd_degree_nodes [2];
            int diff_deg;
            int has_odd_nodes = 0;
            std::shared_ptr<Node> top_node, next_node;

            for (std::unordered_map<std::string, std::shared_ptr<Node>>::iterator it = this->vertices.begin(); it != this->vertices.end(); it++) {
                diff_deg = this->calculate_indegree(it->second) - this->calculate_outdegree(it->second);
                if (diff_deg > 0) {
                    assert(!odd_degree_nodes[0]);
                    odd_degree_nodes[0] = it->second;
                    has_odd_nodes = 1;
                } else if (diff_deg < 0) {
                    assert(!odd_degree_nodes[1]);
                    odd_degree_nodes[1] = it->second;
                    has_odd_nodes = 1;
                }
            }

            assert((odd_degree_nodes[0] && odd_degree_nodes[1]) || (!odd_degree_nodes[0] && !odd_degree_nodes[1]));

            if (has_odd_nodes)
                node_stack.push(odd_degree_nodes[1]);
            else
                node_stack.push(std::next(this->vertices.begin(), rand_between(0, this->vertices.size()))->second);

            while(!node_stack.empty()) {
                top_node = node_stack.top();

                if (!top_node->next.empty()) {
                    next_node = std::next(top_node->next.begin(), rand_between(0, top_node->next.size()))->first;
                    while (has_odd_nodes && next_node == odd_degree_nodes[0] && top_node->next.size() > 1 && this->calculate_indegree(odd_degree_nodes[0]) == 1) {
                        next_node = std::next(top_node->next.begin(), rand_between(0, top_node->next.size()))->first;
                    }

                    this->remove_edge(top_node, next_node);
                    node_stack.push(next_node);
                } else {
                    node_path.push_back(node_stack.top());
                    node_stack.pop();
                }
            }
            return node_path;
        }

        std::vector<std::shared_ptr<Node>> random_walk() {
            std::vector<std::shared_ptr<Node>> node_path;
            std::shared_ptr<Node> top_node, prev_node;
            std::shared_ptr<Node> odd_degree_nodes [2];

            for (std::unordered_map<std::string, std::shared_ptr<Node>>::iterator it = this->vertices.begin(); it != this->vertices.end(); it++) {
                int diff_deg = this->calculate_indegree(it->second) - this->calculate_outdegree(it->second);
                if (diff_deg > 0) {
                    assert(!odd_degree_nodes[0]);
                    odd_degree_nodes[0] = it->second;
                }
                else if (diff_deg < 0) {
                    assert(!odd_degree_nodes[1]);
                    odd_degree_nodes[1] = it->second;
                }
            }

            assert((odd_degree_nodes[0] && odd_degree_nodes[1]) || (!odd_degree_nodes[0] && !odd_degree_nodes[1]));

            if (!odd_degree_nodes[0]) {
                this->add_edge(odd_degree_nodes[0], odd_degree_nodes[1]);
            }

            node_path.push_back(std::next(this->vertices.begin(), rand() % this->vertices.size())->second);
            int num_edges = 0;
            while (!node_path.back()->prev.empty()) {
                top_node = node_path.back();
                prev_node = std::next(top_node->prev.begin(), rand() % top_node->prev.size())->first;
                this->remove_edge(top_node, prev_node);
                node_path.push_back(prev_node);
            }

            return node_path; 
        }

        std::string output_genome(std::vector<std::shared_ptr<Node>> path) {
            std::string genome = "";

            for (std::vector<std::shared_ptr<Node>>::reverse_iterator it = path.rbegin(); it != path.rend(); it++) {
                if (it == path.rend() - 1)
                    genome += (*it)->value;
                else
                    genome += (*it)->value.at(0);
            }
            return genome;
        }

    private:
        std::unordered_map<std::string, std::shared_ptr<Node>> vertices;
};

int main(int argc, char* argv[]) {
    srand(time(NULL));
    int num_trials = 100;
    int num_walks = 0;
    std::string seq;

    int num_successes = 0;
    std::cout << "Testing likelihood of reconstructing original genome with no repeats" << std::endl;
    for (int i = 0; i < num_trials; i++) {
        seq = generate_random_string(1000);
        Graph g(seq, 10);
        
        if (g.output_genome(g.find_eulerian_path()) == seq)
            num_successes += 1;
    }
    std::cout << "Likelihood of reconstructing original genome: " << (float)num_successes / num_trials << std::endl;
    std::cout << "==================================" << std::endl;

    std::cout << "Testing likelihood of generating original genome with 10 repeats of length 20; also including individual number of random walks required to generate genome." << std::endl;
    num_successes = 0;
    for (int i = 0; i < num_trials; i++) {
        seq = add_repeats(generate_random_string(1000), generate_random_string(20), 10);
        
        num_walks = 0;
        while (true) {
            Graph g(seq, 23);
            num_walks += 1;
            if (g.output_genome(g.find_eulerian_path()) == seq) {
                if (num_walks == 1) {
                    num_successes += 1;
                }
                std::cout << num_walks << " ";
                break;
            }
        }
    }
    
    std::cout << std::endl << "Likelihood of reconstructing original genome: " << (float)num_successes / num_trials << std::endl;
    std::cout << "==================================" << std::endl;
    std::cout << "Varying k-mer length" << std::endl;
    num_trials = 10;
    
    int repeat_length = 20;
    int num_repeats = 10;
    for(int read_length = 23; read_length < 30; read_length++) {
        num_walks = 0;
        for (int i = 0; i < num_trials; i++) {
            seq = add_repeats(generate_random_string(1000), generate_random_string(repeat_length), num_repeats);
            int tmp_num_walks = 0;
            while (true) {
                tmp_num_walks++;
                Graph g(seq, read_length);
                if (tmp_num_walks >= 10000 || g.output_genome(g.find_eulerian_path()) == seq) {
                    break;
                }
            }
            num_walks += tmp_num_walks;
        }

        std::cout << "read_length: " << read_length << ", repeat_length: " << repeat_length << ", num_repeats: " << num_repeats << "; num_walks: " << float(num_walks) / num_trials << std::endl;
    }
    std::cout << "==================================" << std::endl;
    std::cout << "Varying repeat number" << std::endl;

    repeat_length = 20;
    int read_length = 10;
    for(int num_repeats = 2; num_repeats < 9; num_repeats++) {
        num_walks = 0;
        for (int i = 0; i < num_trials; i++) {
            seq = add_repeats(generate_random_string(1000), generate_random_string(repeat_length), num_repeats);
            int tmp_num_walks = 0;
            while (true) {
                tmp_num_walks++;
                Graph g(seq, read_length);
                if (tmp_num_walks >= 10000 || g.output_genome(g.find_eulerian_path()) == seq) {
                    break;
                }
            }
            num_walks += tmp_num_walks;
        }

        std::cout << "read_length: " << read_length << ", repeat_length: " << repeat_length << ", num_repeats: " << num_repeats << "; num_walks: " << float(num_walks) / num_trials << std::endl;
    }

    std::cout << "==================================" << std::endl;
    std::cout << "varying repeat length" << std::endl;

    read_length = 10;
    num_repeats = 10;
    for(int repeat_length = 2; repeat_length < 9; repeat_length++) {
        num_walks = 0;
        for (int i = 0; i < num_trials; i++) {
            seq = add_repeats(generate_random_string(1000), generate_random_string(repeat_length), num_repeats);
            int tmp_num_walks = 0;
            while (true) {
                tmp_num_walks++;
                Graph g(seq, read_length);
                if (tmp_num_walks >= 10000 || g.output_genome(g.find_eulerian_path()) == seq) {
                    break;
                }
            }
            num_walks += tmp_num_walks;
        }

        std::cout << "read_length: " << read_length << ", repeat_length: " << repeat_length << ", num_repeats: " << num_repeats << "; num_walks: " << float(num_walks) / num_trials << std::endl;
    }
    
    return 0;
}