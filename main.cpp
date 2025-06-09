#include <bits/stdc++.h>
#include <windows.h>
#include <direct.h>
using namespace std;

const string SCHEMATICS_DIR = "C:\\Users\\Erfan\\Dev\\Cpp\\SUTSpice\\schematics\\";

inline string to_lower_util(string s) {
    transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c){ return tolower(c); });
    return s;
}

string trim_string_util(const string& str) {
    const string whitespace = " \t\n\r\f\v";
    size_t start = str.find_first_not_of(whitespace);
    if (start == string::npos) {
        return "";
    }
    size_t end = str.find_last_not_of(whitespace);
    return str.substr(start, end - start + 1);
}

double parse_value_with_metric_prefix_util(const string& val_str_orig) {
    string val_str = trim_string_util(val_str_orig);
    if (val_str.empty()) {
        throw invalid_argument("Empty value string for parsing.");
    }

    size_t suffix_start_pos = 0;
    double base_val;

    try {
        base_val = stod(val_str, & suffix_start_pos);
    } catch (const invalid_argument & ) {
        throw invalid_argument("Invalid numeric value in string: \"" + val_str_orig + "\"");
    } catch (const out_of_range & ) {
        throw out_of_range("Numeric value out of range in string: \"" + val_str_orig + "\"");
    }

    if (suffix_start_pos >= val_str.length()) {
        return base_val;
    }

    string suffix_str = to_lower_util(val_str.substr(suffix_start_pos));

    double multiplier = 1.0;

    if (suffix_str.rfind("meg", 0) == 0) {
        multiplier = 1e6;
    } else if (!suffix_str.empty()) {
        char prefix_char = suffix_str[0];
        switch (prefix_char) {
            case 'p': multiplier = 1e-12; break;
            case 'n': multiplier = 1e-9;  break;
            case 'u': multiplier = 1e-6;  break;
            case 'm': multiplier = 1e-3;  break;
            case 'k': multiplier = 1e3;   break;
            case 'g': multiplier = 1e9;   break;
        }
    }

    return base_val * multiplier;
}


vector < double > gaussian_elimination_matrix(vector < vector < double >> A, vector < double > b) {
    int n = A.size();
    if (n == 0 || (n > 0 && (A[0].size() != n || b.size() != n))) {
        throw runtime_error("Invalid matrix or vector dimensions for Gaussian elimination.");
    }

    for (int i = 0; i < n; ++i) {
        int max_row = i;
        for (int k = i + 1; k < n; ++k) {
            if (abs(A[k][i]) > abs(A[max_row][i])) max_row = k;
        }

        swap(A[i], A[max_row]);
        swap(b[i], b[max_row]);

        if (abs(A[i][i]) < 1e-12) {
        }

        for (int k = i + 1; k < n; ++k) {
            if (abs(A[i][i]) < 1e-12) continue;
            double factor = A[k][i] / A[i][i];
            for (int j = i; j < n; ++j) A[k][j] -= factor * A[i][j];
            b[k] -= factor * b[i];
        }
    }

    vector < double > x(n);
    for (int i = n - 1; i >= 0; --i) {
        if (abs(A[i][i]) < 1e-12) {
            if (abs(b[i]) > 1e-12) {
                throw runtime_error("System is inconsistent (no solution).");
            }
            x[i] = 0;
        } else {
            x[i] = b[i];
            for (int j = i + 1; j < n; ++j) x[i] -= A[i][j] * x[j];
            x[i] /= A[i][i];
        }
    }
    return x;
}


class Circuit;

enum class ComponentType {
    Resistor,
    Capacitor,
    Inductor,
    VoltageSource,
    CurrentSource,
    VCVS,
    VCCS,
    CCVS,
    CCCS,
    Diode
};
using ResultPoint = map < string, double > ;
const string GROUND_NODE_NAME_CONST = "0";

class Component {
public: string name;
    string node1_name,
            node2_name;
    double value;
    Component(string name_val, string n1, string n2, double val = 0.0): name(move(name_val)),
                                                                                       node1_name(move(n1)),
                                                                                       node2_name(move(n2)),
                                                                                       value(val) {}
    Component(string name_val, string n1, string n2,
              const string & val_str): name(move(name_val)),
                                            node1_name(move(n1)),
                                            node2_name(move(n2)) {
        this -> value = parse_value_with_metric_prefix_util(val_str);
    }
    virtual~Component() =
    default;
    virtual ComponentType get_type() const = 0;
    virtual string to_netlist_string() const = 0;
    virtual void stamp(Circuit & , vector < vector < double >> & , vector < vector < double >> & , vector < vector < double >> & , vector < vector < double >> & , vector < double > & , vector < double > & , map < string, int > & , double,
                       const vector < double > & ) = 0;
    virtual void update_time_dependant_value(double) {}
};

class Resistor: public Component {
public: Resistor(const string & r_name,
                 const string & n1,
                 const string & n2,
                 const string & val_str): Component(r_name, n1, n2, val_str) {
        if (value <= 0) throw runtime_error("Resistor " + name + " must have positive resistance.");
    }
    ComponentType get_type() const override {
        return ComponentType::Resistor;
    }
    string to_netlist_string() const override {
        ostringstream oss;
        oss << name << " " << node1_name << " " << node2_name << " " << fixed << setprecision(12) << value;
        return oss.str();
    }
    void stamp(Circuit & circuit, vector < vector < double >> & G, vector < vector < double >> & , vector < vector < double >> & , vector < vector < double >> & , vector < double > & , vector < double > & , map < string, int > & , double,
               const vector < double > & ) override;
};

class VoltageSource: public Component {
public: enum SourceType {
        DC,
        SINUSOIDAL,
        PULSE
    };
    SourceType sourceType;
    double dc_offset,
            amplitude,
            frequency;
    double v1,
            v2,
            td,
            tr,
            tf,
            pw,
            per;
    vector < string > raw_params;
    VoltageSource(const string & v_name,
                  const string & n1,
                  const string & n2,
                  const vector < string > & params);
    ComponentType get_type() const override {
        return ComponentType::VoltageSource;
    }
    string to_netlist_string() const override;
    void stamp(Circuit & circuit, vector < vector < double >> & , vector < vector < double >> & B, vector < vector < double >> & C, vector < vector < double >> & , vector < double > & , vector < double > & E, map < string, int > & m_map, double,
               const vector < double > & ) override;
    void update_time_dependant_value(double time) override;
};

class CurrentSource: public Component {
public: CurrentSource(const string & i_name,
                      const string & n1,
                      const string & n2,
                      const string & val_str): Component(i_name, n1, n2, val_str) {}
    ComponentType get_type() const override {
        return ComponentType::CurrentSource;
    }
    string to_netlist_string() const override {
        ostringstream oss;
        oss << name << " " << node1_name << " " << node2_name << " " << fixed << setprecision(12) << value;
        return oss.str();
    }
    void stamp(Circuit & circuit, vector < vector < double >> & , vector < vector < double >> & , vector < vector < double >> & , vector < vector < double >> & , vector < double > & J, vector < double > & , map < string, int > & , double,
               const vector < double > & ) override;
};

class Capacitor: public Component {
public: Capacitor(const string & c_name,
                  const string & n1,
                  const string & n2,
                  const string & val_str): Component(c_name, n1, n2, val_str) {
        if (value <= 0) throw runtime_error("Capacitance must be positive.");
    }
    ComponentType get_type() const override {
        return ComponentType::Capacitor;
    }
    string to_netlist_string() const override {
        ostringstream oss;
        oss << name << " " << node1_name << " " << node2_name << " " << fixed << setprecision(12) << value;
        return oss.str();
    }
    void stamp(Circuit & circuit, vector < vector < double >> & G, vector < vector < double >> & , vector < vector < double >> & , vector < vector < double >> & , vector < double > & J, vector < double > & , map < string, int > & , double h,
               const vector < double > & prev_sol) override;
};

class Inductor: public Component {
public: Inductor(const string & l_name,
                 const string & n1,
                 const string & n2,
                 const string & val_str): Component(l_name, n1, n2, val_str) {
        if (value <= 0) throw runtime_error("Inductance must be positive.");
    }
    ComponentType get_type() const override {
        return ComponentType::Inductor;
    }
    string to_netlist_string() const override {
        ostringstream oss;
        oss << name << " " << node1_name << " " << node2_name << " " << fixed << setprecision(12) << value;
        return oss.str();
    }
    void stamp(Circuit & circuit, vector < vector < double >> & , vector < vector < double >> & B, vector < vector < double >> & C, vector < vector < double >> & D, vector < double > & , vector < double > & E, map < string, int > & m_map, double h,
               const vector < double > & prev_sol) override;
};

class Diode: public Component {
public:
    bool is_on = false;
    const double Ron = 1e-3;
    const double Roff = 1e9;

    Diode(const string& name, const string& n1, const string& n2, const string& model)
            : Component(name, n1, n2, 0.0) {
        if (to_lower_util(model) != "ideal") {
            throw runtime_error("Only 'ideal' diode model is supported. Got: " + model);
        }
    }

    ComponentType get_type() const override {
        return ComponentType::Diode;
    }

    string to_netlist_string() const override {
        return name + " " + node1_name + " " + node2_name + " ideal";
    }

    void stamp(Circuit& circuit, vector<vector<double>>& G, vector<vector<double>>& B, vector<vector<double>>& C, vector<vector<double>>& D, vector<double>& J, vector<double>& E, map<string, int>& m_map, double h, const vector<double>& prev_sol) override;
};


class VCVS: public Component {
public: string ctrl_node1,
            ctrl_node2;
    VCVS(const string & name,
         const string & n1,
         const string & n2,
         const string & cn1,
         const string & cn2,
         const string & gain_str): Component(name, n1, n2, gain_str),
                                        ctrl_node1(cn1),
                                        ctrl_node2(cn2) {}
    ComponentType get_type() const override {
        return ComponentType::VCVS;
    }
    string to_netlist_string() const override {
        ostringstream oss;
        oss << name << " " << node1_name << " " << node2_name << " " << ctrl_node1 << " " << ctrl_node2 << " " << fixed << setprecision(12) << value;
        return oss.str();
    }
    void stamp(Circuit & circuit, vector < vector < double >> & , vector < vector < double >> & B, vector < vector < double >> & C, vector < vector < double >> & , vector < double > & , vector < double > & , map < string, int > & m_map, double,
               const vector < double > & ) override;
};

class VCCS: public Component {
public: string ctrl_node1,
            ctrl_node2;
    VCCS(const string & name,
         const string & n1,
         const string & n2,
         const string & cn1,
         const string & cn2,
         const string & gain_str): Component(name, n1, n2, gain_str),
                                        ctrl_node1(cn1),
                                        ctrl_node2(cn2) {}
    ComponentType get_type() const override {
        return ComponentType::VCCS;
    }
    string to_netlist_string() const override {
        ostringstream oss;
        oss << name << " " << node1_name << " " << node2_name << " " << ctrl_node1 << " " << ctrl_node2 << " " << fixed << setprecision(12) << value;
        return oss.str();
    }
    void stamp(Circuit & circuit, vector < vector < double >> & G, vector < vector < double >> & , vector < vector < double >> & , vector < vector < double >> & , vector < double > & , vector < double > & , map < string, int > & , double,
               const vector < double > & ) override;
};

class CCVS: public Component {
public: string ctrl_v_name;
    CCVS(const string & name,
         const string & n1,
         const string & n2,
         const string & cvn,
         const string & gain_str): Component(name, n1, n2, gain_str),
                                        ctrl_v_name(cvn) {}
    ComponentType get_type() const override {
        return ComponentType::CCVS;
    }
    string to_netlist_string() const override {
        ostringstream oss;
        oss << name << " " << node1_name << " " << node2_name << " " << ctrl_v_name << " " << fixed << setprecision(12) << value;
        return oss.str();
    }
    void stamp(Circuit & circuit, vector < vector < double >> & , vector < vector < double >> & B, vector < vector < double >> & C, vector < vector < double >> & D, vector < double > & , vector < double > & , map < string, int > & m_map, double,
               const vector < double > & ) override;
};

class CCCS: public Component {
public: string ctrl_v_name;
    CCCS(const string & name,
         const string & n1,
         const string & n2,
         const string & cvn,
         const string & gain_str): Component(name, n1, n2, gain_str),
                                        ctrl_v_name(cvn) {}
    ComponentType get_type() const override {
        return ComponentType::CCCS;
    }
    string to_netlist_string() const override {
        ostringstream oss;
        oss << name << " " << node1_name << " " << node2_name << " " << ctrl_v_name << " " << fixed << setprecision(12) << value;
        return oss.str();
    }
    void stamp(Circuit & circuit, vector < vector < double >> & , vector < vector < double >> & B, vector < vector < double >> & , vector < vector < double >> & , vector < double > & , vector < double > & , map < string, int > & m_map, double,
               const vector < double > & ) override;
};

class Circuit {
public: vector < unique_ptr < Component >> components;
    map < string,
            int > node_to_idx;
    vector < string > idx_to_node_name;
    string ground_node_explicit_name = GROUND_NODE_NAME_CONST;
    bool ground_node_exists = false;
    vector < VoltageSource * > voltage_source_list;
    vector < Inductor * > inductor_list;
    vector < VCVS * > vcvs_list;
    vector < CCVS * > ccvs_list;
    vector < Diode * > diode_list;
    bool tran_solved = false;
    vector < ResultPoint > tran_results;
    bool dc_sweep_solved = false;
    vector < ResultPoint > dc_sweep_results;
    string dc_sweep_source_name;
    string current_filepath;
    bool is_dirty = false;

    Circuit() =
    default;

    void clear() {
        components.clear();
        node_to_idx.clear();
        idx_to_node_name.clear();
        ground_node_explicit_name = GROUND_NODE_NAME_CONST;
        ground_node_exists = false;
        voltage_source_list.clear();
        inductor_list.clear();
        vcvs_list.clear();
        ccvs_list.clear();
        diode_list.clear();
        tran_solved = false;
        tran_results.clear();
        dc_sweep_solved = false;
        dc_sweep_results.clear();
        dc_sweep_source_name = "";
        current_filepath = "";
        is_dirty = false;
        cout << "INFO: Circuit has been cleared." << endl;
    }

    void save_to_file(const string & path) {
        ofstream file(path);
        if (!file.is_open()) {
            throw runtime_error("Error: Could not open file " + path + " to save.");
        }
        for (const auto & comp: components) {
            file << comp -> to_netlist_string() << "\n";
        }
        if (ground_node_exists) {
            file << "add GND " << ground_node_explicit_name << "\n";
        }
        is_dirty = false;
        current_filepath = path;
        cout << "INFO: Circuit successfully saved to " << path << endl;
    }

    void set_ground_node(const string & node_name) {
        ground_node_explicit_name = node_name;
        ground_node_exists = true;
        is_dirty = true;
        cout << "INFO: Ground node set to/confirmed as: " << ground_node_explicit_name << endl;
    }

    void delete_ground_node(const string & node_name) {
        if (!ground_node_exists || node_name != ground_node_explicit_name) throw runtime_error("Node does not exist");
        ground_node_exists = false;
        ground_node_explicit_name = GROUND_NODE_NAME_CONST;
        tran_solved = false;
        dc_sweep_solved = false;
        is_dirty = true;
        cout << "Ground reference on node '" << node_name << "' removed." << endl;
    }

    bool is_ground(const string & node_name) const {
        if (!ground_node_exists) return false;
        if (node_name == ground_node_explicit_name) return true;
        if (ground_node_explicit_name == "0" && node_name == "GND") return true;
        if (ground_node_explicit_name == "GND" && node_name == "0") return true;
        return false;
    }

    void add_component(unique_ptr < Component > comp) {
        for (const auto & existing_comp: components) {
            if (existing_comp -> name == comp -> name) throw runtime_error("Error: Component with name '" + comp -> name + "' already exists.");
        }
        tran_solved = false;
        dc_sweep_solved = false;
        is_dirty = true;
        if (comp -> node1_name == GROUND_NODE_NAME_CONST || comp -> node1_name == "GND" ||
            comp -> node2_name == GROUND_NODE_NAME_CONST || comp -> node2_name == "GND") {
            if (!ground_node_exists) {
                ground_node_explicit_name = (comp -> node1_name == "GND" || comp -> node2_name == "GND") ? "GND" : GROUND_NODE_NAME_CONST;
                ground_node_exists = true;
            }
        }
        components.push_back(move(comp));
    }

    void delete_component(const string & component_name) {
        auto it = find_if(components.begin(), components.end(), [ & ](const auto & comp) {
            return comp -> name == component_name;
        });
        if (it == components.end()) throw runtime_error("Cannot delete component; component '" + component_name + "' not found.");
        components.erase(it);
        tran_solved = false;
        dc_sweep_solved = false;
        is_dirty = true;
        cout << "Component '" << component_name << "' deleted successfully." << endl;
    }

    void rename_node(const string & old_name,
                     const string & new_name) {
        set < string > all_nodes;
        for (const auto & comp: components) {
            all_nodes.insert(comp -> node1_name);
            all_nodes.insert(comp -> node2_name);
        }
        if (all_nodes.find(old_name) == all_nodes.end()) throw runtime_error("ERROR: Node " + old_name + " does not exist in the circuit");
        if (all_nodes.find(new_name) != all_nodes.end()) throw runtime_error("ERROR: Node name " + new_name + " already exists");
        for (auto & comp: components) {
            if (comp -> node1_name == old_name) comp -> node1_name = new_name;
            if (comp -> node2_name == old_name) comp -> node2_name = new_name;
        }
        if (ground_node_exists && ground_node_explicit_name == old_name) ground_node_explicit_name = new_name;
        tran_solved = false;
        dc_sweep_solved = false;
        is_dirty = true;
        cout << "SUCCESS: Node renamed from " << old_name << " to " << new_name << endl;
    }

    int prepare_for_analysis() {
        node_to_idx.clear();
        idx_to_node_name.clear();
        voltage_source_list.clear();
        inductor_list.clear();
        vcvs_list.clear();
        ccvs_list.clear();
        diode_list.clear();
        set < string > unique_node_names;
        for (const auto & comp: components) {
            unique_node_names.insert(comp -> node1_name);
            unique_node_names.insert(comp -> node2_name);
            if (auto vs = dynamic_cast < VoltageSource * > (comp.get())) voltage_source_list.push_back(vs);
            else if (auto ind = dynamic_cast < Inductor * > (comp.get())) inductor_list.push_back(ind);
            else if (auto vcvs = dynamic_cast < VCVS * > (comp.get())) vcvs_list.push_back(vcvs);
            else if (auto ccvs = dynamic_cast < CCVS * > (comp.get())) ccvs_list.push_back(ccvs);
            else if (auto d = dynamic_cast < Diode * > (comp.get())) diode_list.push_back(d);
        }
        if (!ground_node_exists && (unique_node_names.count("0") || unique_node_names.count("GND"))) {
            set_ground_node(unique_node_names.count("0") ? "0" : "GND");
        }
        int idx = 0;
        for (const auto & name: unique_node_names)
            if (!is_ground(name)) {
                node_to_idx[name] = idx++;
                idx_to_node_name.push_back(name);
            }
        return idx;
    }

    int get_node_matrix_index(const string & node_name) const {
        if (is_ground(node_name)) return -1;
        auto it = node_to_idx.find(node_name);
        if (it != node_to_idx.end()) {
            return it -> second;
        }
        return -2;
    }

    double get_voltage_at(const string & node_name,
                          const ResultPoint & results) const {
        if (is_ground(node_name)) return 0.0;
        string v_name = "V(" + node_name + ")";
        auto it = results.find(v_name);
        if (it != results.end()) return it -> second;
        return 0.0;
    }

    void calculate_and_store_passive_currents(ResultPoint & result_point,
                                              const ResultPoint & prev_result_point, double h) {
        for (const auto & comp: components) {
            double v1 = get_voltage_at(comp -> node1_name, result_point);
            double v2 = get_voltage_at(comp -> node2_name, result_point);

            if (auto r = dynamic_cast < Resistor * > (comp.get())) {
                result_point["I(" + r -> name + ")"] = (v1 - v2) / r -> value;
            } else if (auto c = dynamic_cast < Capacitor * > (comp.get())) {
                if (h > 0 && !prev_result_point.empty()) {
                    double v1_prev = get_voltage_at(c -> node1_name, prev_result_point);
                    double v2_prev = get_voltage_at(c -> node2_name, prev_result_point);
                    result_point["I(" + c -> name + ")"] = c -> value * ((v1 - v2) - (v1_prev - v2_prev)) / h;
                } else {
                    result_point["I(" + c -> name + ")"] = 0.0;
                }
            } else if (auto d = dynamic_cast<Diode*>(comp.get())) {
                result_point["I(" + d->name + ")"] = (v1 - v2) / (d->is_on ? d->Ron : d->Roff);
            }
        }
    }

    void build_mna_matrix(vector < vector < double >> & A, vector < double > & z, double h,
                          const vector < double > & prev_sol) {
        if (!components.empty() && !ground_node_exists) throw runtime_error("No ground node defined.");
        int N = idx_to_node_name.size();
        int M = voltage_source_list.size() + inductor_list.size() + vcvs_list.size() + ccvs_list.size();
        int system_size = N + M;
        if (system_size == 0 && diode_list.empty()) {
            A.clear();
            z.clear();
            return;
        }
        A.assign(system_size, vector < double > (system_size, 0.0));
        z.assign(system_size, 0.0);
        vector < vector < double >> G(N, vector < double > (N, 0.0)), B(N, vector < double > (M, 0.0)), C(M, vector < double > (N, 0.0)), D(M, vector < double > (M, 0.0));
        vector < double > J(N, 0.0), E(M, 0.0);
        map < string, int > m_map;
        int m_counter = 0;
        for (const auto & vs: voltage_source_list) m_map[vs -> name] = m_counter++;
        for (const auto & l: inductor_list) m_map[l -> name] = m_counter++;
        for (const auto & e: vcvs_list) m_map[e -> name] = m_counter++;
        for (const auto & h: ccvs_list) m_map[h -> name] = m_counter++;
        for (const auto & comp: components) comp -> stamp( * this, G, B, C, D, J, E, m_map, h, prev_sol);
        for (int i = 0; i < N; ++i)
            for (int j = 0; j < N; ++j) A[i][j] = G[i][j];
        for (int i = 0; i < N; ++i)
            for (int j = 0; j < M; ++j) A[i][N + j] = B[i][j];
        for (int i = 0; i < M; ++i)
            for (int j = 0; j < N; ++j) A[N + i][j] = C[i][j];
        for (int i = 0; i < M; ++i)
            for (int j = 0; j < M; ++j) A[N + i][N + j] = D[i][j];
        for (int i = 0; i < N; ++i) z[i] = J[i];
        for (int i = 0; i < M; ++i) z[N + i] = E[i];
    }

    void perform_transient_analysis(double t_step, double t_stop) {
        tran_solved = false;
        tran_results.clear();
        if (t_step <= 0 || t_stop <= 0 || t_step > t_stop) throw runtime_error("Invalid transient parameters.");

        int N = prepare_for_analysis();
        int M = voltage_source_list.size() + inductor_list.size() + vcvs_list.size() + ccvs_list.size();

        if (N + M == 0 && diode_list.empty()) {
            tran_solved = true;
            return;
        }

        for (auto* d : diode_list) {
            d->is_on = false;
        }

        vector<double> prev_solution(N + M, 0.0);
        ResultPoint prev_result_point;

        for (double t = 0; t <= t_stop + (t_step / 2.0); t += t_step) {
            for (auto& comp : components) {
                comp->update_time_dependant_value(t);
            }

            const int MAX_DIODE_ITERATIONS = 100;
            bool diodes_converged = false;
            vector<double> current_solution;

            for (int iter = 0; iter < MAX_DIODE_ITERATIONS; ++iter) {
                vector<vector<double>> A;
                vector<double> z;
                build_mna_matrix(A, z, t_step, prev_solution);
                current_solution = gaussian_elimination_matrix(A, z);

                if (diode_list.empty()) {
                    diodes_converged = true;
                    break;
                }

                bool state_changed = false;
                for (auto* d : diode_list) {
                    bool old_state = d->is_on;
                    double v1 = (get_node_matrix_index(d->node1_name) >= 0) ? current_solution[get_node_matrix_index(d->node1_name)] : 0.0;
                    double v2 = (get_node_matrix_index(d->node2_name) >= 0) ? current_solution[get_node_matrix_index(d->node2_name)] : 0.0;

                    d->is_on = (v1 > v2);

                    if (d->is_on != old_state) {
                        state_changed = true;
                    }
                }

                if (!state_changed) {
                    diodes_converged = true;
                    break;
                }
            }

            if (!diodes_converged) {
                throw runtime_error("Diode model failed to converge at time " + to_string(t));
            }

            ResultPoint result_at_t;
            result_at_t["time"] = t;
            for (int i = 0; i < N; ++i) result_at_t["V(" + idx_to_node_name[i] + ")"] = current_solution[i];
            if (ground_node_exists) result_at_t["V(" + ground_node_explicit_name + ")"] = 0.0;

            map<string, int> m_map;
            int m_counter = 0;
            for (const auto & vs: voltage_source_list) m_map[vs -> name] = m_counter++;
            for (const auto & l: inductor_list) m_map[l -> name] = m_counter++;
            for (const auto & e: vcvs_list) m_map[e -> name] = m_counter++;
            for (const auto & h: ccvs_list) m_map[h -> name] = m_counter++;

            for (const auto & p: m_map) result_at_t["I(" + p.first + ")"] = current_solution[N + p.second];

            calculate_and_store_passive_currents(result_at_t, prev_result_point, t_step);
            tran_results.push_back(result_at_t);
            prev_solution = current_solution;
            prev_result_point = result_at_t;
        }
        tran_solved = true;
        cout << "Transient analysis completed." << endl;
    }

    void perform_dc_sweep_analysis(const string & src, double start, double end, double inc) {
        dc_sweep_solved = false;
        dc_sweep_results.clear();
        if (inc == 0 || (end > start && inc < 0) || (end < start && inc > 0)) throw runtime_error("Invalid sweep parameters.");
        Component * sweep_comp = nullptr;
        for (auto & c: components)
            if (c -> name == src) {
                sweep_comp = c.get();
                break;
            }
        if (!sweep_comp || (sweep_comp -> get_type() != ComponentType::VoltageSource && sweep_comp -> get_type() != ComponentType::CurrentSource)) throw runtime_error("Sweep source not found or invalid.");

        double orig_val = sweep_comp -> value;
        dc_sweep_source_name = src;
        int N = prepare_for_analysis();
        int M = voltage_source_list.size() + inductor_list.size() + vcvs_list.size() + ccvs_list.size();
        if (N + M == 0 && diode_list.empty()) {
            dc_sweep_solved = true;
            return;
        }

        for (auto* d : diode_list) {
            d->is_on = false;
        }

        for (double val = start;
             (inc > 0) ? (val <= end + abs(inc) / 2.0) : (val >= end - abs(inc) / 2.0); val += inc) {
            sweep_comp -> value = val;

            const int MAX_DIODE_ITERATIONS = 100;
            bool diodes_converged = false;
            vector<double> solution;

            for (int iter = 0; iter < MAX_DIODE_ITERATIONS; ++iter) {
                vector<vector<double>> A;
                vector<double> z;
                build_mna_matrix(A, z, 0.0, {});
                solution = gaussian_elimination_matrix(A, z);

                if (diode_list.empty()) {
                    diodes_converged = true;
                    break;
                }

                bool state_changed = false;
                for (auto* d : diode_list) {
                    bool old_state = d->is_on;
                    double v1 = (get_node_matrix_index(d->node1_name) >= 0 && get_node_matrix_index(d->node1_name) < solution.size()) ? solution[get_node_matrix_index(d->node1_name)] : 0.0;
                    double v2 = (get_node_matrix_index(d->node2_name) >= 0 && get_node_matrix_index(d->node2_name) < solution.size()) ? solution[get_node_matrix_index(d->node2_name)] : 0.0;

                    d->is_on = (v1 > v2);

                    if (d->is_on != old_state) {
                        state_changed = true;
                    }
                }

                if (!state_changed) {
                    diodes_converged = true;
                    break;
                }
            }

            if (!diodes_converged) {
                throw runtime_error("Diode model failed to converge at sweep value " + to_string(val));
            }

            ResultPoint result_at_val;
            result_at_val[src] = val;
            for (int i = 0; i < N; ++i) result_at_val["V(" + idx_to_node_name[i] + ")"] = solution[i];
            if (ground_node_exists) result_at_val["V(" + ground_node_explicit_name + ")"] = 0.0;

            map < string, int > m_map;
            int m_counter = 0;
            for (const auto & vs: voltage_source_list) m_map[vs -> name] = m_counter++;
            for (const auto & l: inductor_list) m_map[l -> name] = m_counter++;
            for (const auto & e: vcvs_list) m_map[e -> name] = m_counter++;
            for (const auto & h: ccvs_list) m_map[h -> name] = m_counter++;

            for (const auto & p: m_map) result_at_val["I(" + p.first + ")"] = solution[N + p.second];

            calculate_and_store_passive_currents(result_at_val, {}, 0.0);
            dc_sweep_results.push_back(result_at_val);
        }
        sweep_comp -> value = orig_val;
        dc_sweep_solved = true;
        cout << "DC sweep analysis completed." << endl;
    }
};

VoltageSource::VoltageSource(const string & v_name,
                             const string & n1,
                             const string & n2,
                             const vector < string > & params): Component(v_name, n1, n2, 0.0), sourceType(DC), dc_offset(0.0), amplitude(0.0), frequency(0.0), v1(0), v2(0), td(0), tr(0), tf(0), pw(0), per(0), raw_params(params) {
    if (params.empty()) throw runtime_error("No value or parameters provided for voltage source " + name);
    string first_param_upper = to_lower_util(params[0]);
    if (first_param_upper.rfind("sin", 0) == 0) {
        string combined_params;
        for (size_t i = 0; i < params.size(); ++i) {
            combined_params += params[i];
            if (i < params.size() - 1) combined_params += " ";
        }
        size_t start_paren = combined_params.find('('), end_paren = combined_params.rfind(')');
        if (start_paren == string::npos || end_paren == string::npos) throw runtime_error("Mismatched parentheses in SIN() for " + name);
        string sin_params_str = combined_params.substr(start_paren + 1, end_paren - start_paren - 1);
        stringstream ss(sin_params_str);
        string offset_str, amp_str, freq_str;
        ss >> offset_str >> amp_str >> freq_str;
        if (ss.fail() || !ss.eof()) throw runtime_error("Invalid parameters inside SIN() for " + name);
        sourceType = SINUSOIDAL;
        dc_offset = parse_value_with_metric_prefix_util(offset_str);
        amplitude = parse_value_with_metric_prefix_util(amp_str);
        frequency = parse_value_with_metric_prefix_util(freq_str);
        this -> value = dc_offset;
    } else if (first_param_upper.rfind("pulse", 0) == 0) {
        string combined_params;
        for (size_t i = 0; i < params.size(); ++i) {
            combined_params += params[i];
            if (i < params.size() - 1) combined_params += " ";
        }
        size_t start_paren = combined_params.find('('), end_paren = combined_params.rfind(')');
        if (start_paren == string::npos || end_paren == string::npos) throw runtime_error("Mismatched parentheses in PULSE() for " + name);
        string pulse_params_str = combined_params.substr(start_paren + 1, end_paren - start_paren - 1);
        stringstream ss(pulse_params_str);
        string v1_str, v2_str, td_str, tr_str, tf_str, pw_str, per_str;
        ss >> v1_str >> v2_str >> td_str >> tr_str >> tf_str >> pw_str >> per_str;
        if (ss.fail() || !ss.eof()) throw runtime_error("Invalid parameters inside PULSE() for " + name);
        sourceType = PULSE;
        v1 = parse_value_with_metric_prefix_util(v1_str);
        v2 = parse_value_with_metric_prefix_util(v2_str);
        td = parse_value_with_metric_prefix_util(td_str);
        tr = parse_value_with_metric_prefix_util(tr_str);
        tf = parse_value_with_metric_prefix_util(tf_str);
        pw = parse_value_with_metric_prefix_util(pw_str);
        per = parse_value_with_metric_prefix_util(per_str);
        this -> value = v1;
    } else if (params.size() == 1) {
        sourceType = DC;
        this -> value = parse_value_with_metric_prefix_util(params[0]);
        this -> dc_offset = this -> value;
    } else throw runtime_error("Invalid parameters for voltage source " + name);
}

string VoltageSource::to_netlist_string() const {
    ostringstream oss;
    oss << name << " " << node1_name << " " << node2_name << " ";
    for (size_t i = 0; i < raw_params.size(); ++i) {
        oss << raw_params[i] << (i == raw_params.size() - 1 ? "" : " ");
    }
    return oss.str();
}

void Resistor::stamp(Circuit & circuit, vector < vector < double >> & G, vector < vector < double >> & B, vector < vector < double >> & C, vector < vector < double >> & D, vector < double > & J, vector < double > & E, map < string, int > & m_map, double h,
                     const vector < double > & prev_sol) {
    double conductance = 1.0 / value;
    int idx1 = circuit.get_node_matrix_index(node1_name);
    int idx2 = circuit.get_node_matrix_index(node2_name);
    if (idx1 >= 0) G[idx1][idx1] += conductance;
    if (idx2 >= 0) G[idx2][idx2] += conductance;
    if (idx1 >= 0 && idx2 >= 0) {
        G[idx1][idx2] -= conductance;
        G[idx2][idx1] -= conductance;
    }
}

void VoltageSource::stamp(Circuit & circuit, vector < vector < double >> & G, vector < vector < double >> & B, vector < vector < double >> & C, vector < vector < double >> & D, vector < double > & J, vector < double > & E, map < string, int > & m_map, double h,
                          const vector < double > & prev_sol) {
    int m_idx = m_map.at(name);
    int idx1 = circuit.get_node_matrix_index(node1_name);
    int idx2 = circuit.get_node_matrix_index(node2_name);
    if (idx1 >= 0) {
        B[idx1][m_idx] += 1.0;
        C[m_idx][idx1] += 1.0;
    }
    if (idx2 >= 0) {
        B[idx2][m_idx] -= 1.0;
        C[m_idx][idx2] -= 1.0;
    }
    E[m_idx] += this -> value;
}

void VoltageSource::update_time_dependant_value(double time) {
    if (sourceType == SINUSOIDAL) {
        if (frequency > 0) {
            value = dc_offset + amplitude * sin(2 * 3.14159265358979323846 * frequency * time);
        } else {
            value = dc_offset + amplitude;
        }
    } else if (sourceType == PULSE) {
        if (time < td) {
            value = v1;
        } else {
            double t_rel = per > 0 ? fmod(time - td, per) : time - td;
            if (tr > 0 && t_rel < tr) value = v1 + (v2 - v1) * (t_rel / tr);
            else if (t_rel < tr + pw) value = v2;
            else if (tf > 0 && t_rel < tr + pw + tf) value = v2 + (v1 - v2) * ((t_rel - tr - pw) / tf);
            else value = v1;
        }
    }
}

void CurrentSource::stamp(Circuit & circuit, vector < vector < double >> & G, vector < vector < double >> & B, vector < vector < double >> & C, vector < vector < double >> & D, vector < double > & J, vector < double > & E, map < string, int > & m_map, double h,
                          const vector < double > & prev_sol) {
    int idx1 = circuit.get_node_matrix_index(node1_name);
    int idx2 = circuit.get_node_matrix_index(node2_name);
    if (idx1 >= 0) J[idx1] -= value;
    if (idx2 >= 0) J[idx2] += value;
}

void Capacitor::stamp(Circuit & circuit, vector < vector < double >> & G, vector < vector < double >> & B, vector < vector < double >> & C, vector < vector < double >> & D, vector < double > & J, vector < double > & E, map < string, int > & m_map, double h,
                      const vector < double > & prev_sol) {
    if (h == 0) return;
    double g_eq = value / h;
    int idx1 = circuit.get_node_matrix_index(node1_name);
    int idx2 = circuit.get_node_matrix_index(node2_name);
    if (idx1 >= 0) G[idx1][idx1] += g_eq;
    if (idx2 >= 0) G[idx2][idx2] += g_eq;
    if (idx1 >= 0 && idx2 >= 0) {
        G[idx1][idx2] -= g_eq;
        G[idx2][idx1] -= g_eq;
    }
    if (!prev_sol.empty()) {
        double v1_prev = (idx1 >= 0) ? prev_sol[idx1] : 0.0;
        double v2_prev = (idx2 >= 0) ? prev_sol[idx2] : 0.0;
        double i_eq = g_eq * (v1_prev - v2_prev);
        if (idx1 >= 0) J[idx1] += i_eq;
        if (idx2 >= 0) J[idx2] -= i_eq;
    }
}

void Inductor::stamp(Circuit & circuit, vector < vector < double >> & G, vector < vector < double >> & B, vector < vector < double >> & C, vector < vector < double >> & D, vector < double > & J, vector < double > & E, map < string, int > & m_map, double h,
                     const vector < double > & prev_sol) {
    int m_idx = m_map.at(name);
    int idx1 = circuit.get_node_matrix_index(node1_name);
    int idx2 = circuit.get_node_matrix_index(node2_name);
    if (idx1 >= 0) {
        B[idx1][m_idx] += 1.0;
        C[m_idx][idx1] += 1.0;
    }
    if (idx2 >= 0) {
        B[idx2][m_idx] -= 1.0;
        C[m_idx][idx2] -= 1.0;
    }
    if (h > 0) {
        double r_eq = value / h;
        D[m_idx][m_idx] -= r_eq;
        int N = circuit.idx_to_node_name.size();
        if (!prev_sol.empty() && (N + m_idx < prev_sol.size())) {
            E[m_idx] -= r_eq * prev_sol[N + m_idx];
        }
    }
}

void Diode::stamp(Circuit& circuit, vector<vector<double>>& G, vector<vector<double>>& B, vector<vector<double>>& C, vector<vector<double>>& D, vector<double>& J, vector<double>& E, map<string, int>& m_map, double h, const vector<double>& prev_sol) {
    double conductance = 1.0 / (is_on ? Ron : Roff);
    int idx1 = circuit.get_node_matrix_index(node1_name);
    int idx2 = circuit.get_node_matrix_index(node2_name);
    if (idx1 >= 0) G[idx1][idx1] += conductance;
    if (idx2 >= 0) G[idx2][idx2] += conductance;
    if (idx1 >= 0 && idx2 >= 0) {
        G[idx1][idx2] -= conductance;
        G[idx2][idx1] -= conductance;
    }
}

void VCVS::stamp(Circuit & circuit, vector < vector < double >> & G, vector < vector < double >> & B, vector < vector < double >> & C, vector < vector < double >> & D, vector < double > & J, vector < double > & E, map < string, int > & m_map, double h,
                 const vector < double > & prev_sol) {
    int m_idx = m_map.at(name);
    int idx_n1 = circuit.get_node_matrix_index(node1_name);
    int idx_n2 = circuit.get_node_matrix_index(node2_name);
    int idx_cn1 = circuit.get_node_matrix_index(ctrl_node1);
    int idx_cn2 = circuit.get_node_matrix_index(ctrl_node2);
    if (idx_n1 >= 0) B[idx_n1][m_idx] = 1.0;
    if (idx_n2 >= 0) B[idx_n2][m_idx] = -1.0;
    if (idx_n1 >= 0) C[m_idx][idx_n1] = 1.0;
    if (idx_n2 >= 0) C[m_idx][idx_n2] = -1.0;
    if (idx_cn1 >= 0) C[m_idx][idx_cn1] -= value;
    if (idx_cn2 >= 0) C[m_idx][idx_cn2] += value;
}

void VCCS::stamp(Circuit & circuit, vector < vector < double >> & G, vector < vector < double >> & B, vector < vector < double >> & C, vector < vector < double >> & D, vector < double > & J, vector < double > & E, map < string, int > & m_map, double h,
                 const vector < double > & prev_sol) {
    int idx_n1 = circuit.get_node_matrix_index(node1_name);
    int idx_n2 = circuit.get_node_matrix_index(node2_name);
    int idx_cn1 = circuit.get_node_matrix_index(ctrl_node1);
    int idx_cn2 = circuit.get_node_matrix_index(ctrl_node2);
    if (idx_n1 >= 0 && idx_cn1 >= 0) G[idx_n1][idx_cn1] += value;
    if (idx_n1 >= 0 && idx_cn2 >= 0) G[idx_n1][idx_cn2] -= value;
    if (idx_n2 >= 0 && idx_cn1 >= 0) G[idx_n2][idx_cn1] -= value;
    if (idx_n2 >= 0 && idx_cn2 >= 0) G[idx_n2][idx_cn2] += value;
}

void CCVS::stamp(Circuit & circuit, vector < vector < double >> & G, vector < vector < double >> & B, vector < vector < double >> & C, vector < vector < double >> & D, vector < double > & J, vector < double > & E, map < string, int > & m_map, double h,
                 const vector < double > & prev_sol) {
    int m_idx_h = m_map.at(name);
    if (m_map.find(ctrl_v_name) == m_map.end()) throw runtime_error("Control source " + ctrl_v_name + " not found for " + name);
    int m_idx_ctrl = m_map.at(ctrl_v_name);
    int idx_n1 = circuit.get_node_matrix_index(node1_name);
    int idx_n2 = circuit.get_node_matrix_index(node2_name);
    if (idx_n1 >= 0) B[idx_n1][m_idx_h] = 1.0;
    if (idx_n2 >= 0) B[idx_n2][m_idx_h] = -1.0;
    if (idx_n1 >= 0) C[m_idx_h][idx_n1] = 1.0;
    if (idx_n2 >= 0) C[m_idx_h][idx_n2] = -1.0;
    D[m_idx_h][m_idx_ctrl] = -value;
}

void CCCS::stamp(Circuit & circuit, vector < vector < double >> & G, vector < vector < double >> & B, vector < vector < double >> & C, vector < vector < double >> & D, vector < double > & J, vector < double > & E, map < string, int > & m_map, double h,
                 const vector < double > & prev_sol) {
    if (m_map.find(ctrl_v_name) == m_map.end()) throw runtime_error("Control source " + ctrl_v_name + " not found for " + name);
    int m_idx_ctrl = m_map.at(ctrl_v_name);
    int idx_n1 = circuit.get_node_matrix_index(node1_name);
    int idx_n2 = circuit.get_node_matrix_index(node2_name);
    if (idx_n1 >= 0) B[idx_n1][m_idx_ctrl] += value;
    if (idx_n2 >= 0) B[idx_n2][m_idx_ctrl] -= value;
}

struct Command {
    string type;
    vector < string > args;
};

class Parser {
public: Command parse_line(const string & line) {
        stringstream ss(trim_string_util(line));
        string segment;
        Command cmd;
        bool first = true;
        while (ss >> segment) {
            if (first) {
                cmd.type = segment;
                first = false;
            } else cmd.args.push_back(segment);
        }
        return cmd;
    }

    void execute_command(Command cmd, Circuit & circuit) {
        try {
            bool circuit_modified = false;
            string cmd_type_lower = to_lower_util(cmd.type);

            set < string > known_commands = {
                    ".nodes",
                    ".list",
                    ".rename",
                    ".print",
                    ".clear",
                    ".show",
                    ".save",
                    ".load",
                    "add",
                    "delete",
                    "exit",
            };
            if (known_commands.find(cmd_type_lower) == known_commands.end() && !cmd.type.empty() && cmd.type[0] != '.') {
                cmd.args.insert(cmd.args.begin(), cmd.type);
                cmd.type = "add";
                cmd_type_lower = "add";
            }

            if (cmd_type_lower == "add") {
                if (cmd.args.empty()) throw runtime_error("Usage: add <item> ...");
                if (cmd.args[0] == "GND") handle_add_ground(cmd.args, circuit);
                else handle_add_component(cmd.args, circuit);
                circuit_modified = true;
            } else if (cmd_type_lower == "delete") {
                if (cmd.args.empty()) throw runtime_error("Usage: delete <item> ...");
                if (cmd.args[0] == "GND") handle_delete_ground(cmd.args, circuit);
                else handle_delete_component(cmd.args, circuit);
                circuit_modified = true;
            } else if (cmd_type_lower == ".nodes") handle_list_nodes(circuit);
            else if (cmd_type_lower == ".list") handle_list_components(cmd.args, circuit);
            else if (cmd_type_lower == ".rename") {
                handle_rename_node(cmd.args, circuit);
                circuit_modified = true;
            } else if (cmd_type_lower == ".print") handle_print(cmd.args, circuit);
            else if (cmd_type_lower == ".clear") circuit.clear();
            else if (cmd_type_lower == ".show") handle_show_schematics(cmd.args, circuit);
            else if (cmd_type_lower == ".save") {
                handle_save(cmd.args, circuit);
            } else if (cmd_type_lower == ".load") handle_load_command(cmd.args, circuit);
            else if (cmd_type_lower == "exit" || cmd.type.empty()) {} else cerr << "Error: Unknown command '" << cmd.type << "'" << endl;

            if (circuit_modified && !circuit.current_filepath.empty()) {
                cout << "INFO: Auto-saving changes to " << circuit.current_filepath << "..." << endl;
                circuit.save_to_file(circuit.current_filepath);
            }

        } catch (const exception & e) {
            cerr << "Error: " << e.what() << endl;
        }
    }
private:
    void handle_load_command(const vector < string > & args, Circuit & circuit) {
        if (args.empty()) {
            cout << "Enter filepath to load: ";
            string path;
            getline(cin, path);
            path = trim_string_util(path);
            if (path.empty()) {
                cout << "Load cancelled." << endl;
                return;
            }
            handle_load_file_from_path({
                                               path
                                       }, circuit);
        } else {
            handle_load_file_from_path(args, circuit);
        }
    }

    void handle_save(const vector<string>& args, Circuit& circuit) {
        string filename_or_path;

        if (args.empty()) {
            cout << "Enter filename to save (or full path): ";
            getline(cin, filename_or_path);
            filename_or_path = trim_string_util(filename_or_path);
            if (filename_or_path.empty()) {
                cout << "Save cancelled." << endl;
                return;
            }
        } else {
            filename_or_path = args[0];
            for (size_t i = 1; i < args.size(); ++i) {
                filename_or_path += " " + args[i];
            }
        }

        string final_path = filename_or_path;
        if (filename_or_path.find('\\') == string::npos && filename_or_path.find('/') == string::npos) {
            final_path = SCHEMATICS_DIR + filename_or_path;
        }

        if (final_path.length() < 4 || to_lower_util(final_path.substr(final_path.length() - 4)) != ".txt") {
            final_path += ".txt";
        }

        circuit.save_to_file(final_path);
    }

    void handle_show_schematics(const vector < string > & args, Circuit & circuit) {
        if (args.size() != 2 || to_lower_util(args[0]) != "existing" || to_lower_util(args[1]) != "schematics") {
            throw runtime_error("Unknown command. Did you mean '.show existing schematics'?");
        }
        while (true) {
            vector < string > schematics;
#ifdef _WIN32
            string search_path = SCHEMATICS_DIR + "*.txt"; WIN32_FIND_DATA fd; HANDLE hFind = FindFirstFile(search_path.c_str(), &fd);
            if (hFind != INVALID_HANDLE_VALUE) { do { schematics.push_back(fd.cFileName); } while (FindNextFile(hFind, &fd) != 0); FindClose(hFind); }
#else
            DIR *dir; struct dirent *ent;
                if ((dir = opendir(SCHEMATICS_DIR.c_str())) != NULL) {
                    while ((ent = readdir(dir)) != NULL) {
                        string fn = ent->d_name;
                        if (fn.length() > 4 && fn.substr(fn.length() - 4) == ".txt") schematics.push_back(fn);
                    }
                    closedir(dir);
                } else throw runtime_error("Could not access directory: " + SCHEMATICS_DIR);
#endif
            cout << "-choose existing schematic:" << endl;
            if (schematics.empty()) cout << "  (No .txt files found in '" << SCHEMATICS_DIR << "' directory)" << endl;
            else
                for (size_t i = 0; i < schematics.size(); ++i) cout << i + 1 << "-" << schematics[i] << endl;
            string choice;
            if (!getline(cin, choice)) break;
            choice = trim_string_util(to_lower_util(choice));
            if (choice == "return") {
                cout << "Returning to main menu." << endl;
                break;
            }
            try {
                int index = stoi(choice);
                if (index < 1 || index > schematics.size()) cout << "-Error : Inappropriate input" << endl;
                else {
                    string filename = schematics[index - 1];
                    handle_load_file_from_path({
                                                       filename
                                               }, circuit);
                    cout << "--- Loaded Circuit Status ---" << endl;
                    handle_list_components({}, circuit);
                    handle_list_nodes(circuit);
                    cout << "---------------------------" << endl;
                    break;
                }
            } catch (const exception & ) {
                cout << "-Error : Inappropriate input" << endl;
            }
        }
    }

    void handle_load_file_from_path(const vector < string > & args, Circuit & circuit) {
        if (args.size() != 1) throw runtime_error("Usage: requires a filename.");
        string filepath = args[0];
        if (filepath.find('/') == string::npos && filepath.find('\\') == string::npos) filepath = SCHEMATICS_DIR + filepath;
        ifstream file(filepath);
        if (!file.is_open()) throw runtime_error("Cannot open file: " + filepath);
        circuit.clear();
        circuit.current_filepath = filepath;
        string line;
        int line_num = 0;
        cout << "INFO: Loading netlist from: " << filepath << endl;
        while (getline(file, line)) {
            line_num++;
            line = trim_string_util(line);
            if (line.empty() || line[0] == '*') continue;
            if (to_lower_util(line) == ".end") {
                cout << "INFO: '.end' directive found. Finished loading." << endl;
                break;
            }
            cout << "  > Executing from file: " << line << endl;
            try {
                Command cmd = parse_line(line);
                string temp_path = circuit.current_filepath;
                circuit.current_filepath = "";
                execute_command(cmd, circuit);
                circuit.current_filepath = temp_path;
            } catch (const exception & e) {
                cerr << "ERROR on line " << line_num << " ('" << line << "'): " << e.what() << endl;
                circuit.clear();
                throw runtime_error("Failed to load file due to error.");
            }
        }
        circuit.is_dirty = false;
    }

    void handle_add_ground(const vector < string > & args, Circuit & circuit) {
        if (args.size() != 2) throw runtime_error("Usage: add GND <node>");
        circuit.set_ground_node(args[1]);
    }
    void handle_delete_ground(const vector < string > & args, Circuit & circuit) {
        if (args.size() != 2) throw runtime_error("Usage: delete GND <node>");
        circuit.delete_ground_node(args[1]);
    }
    void handle_delete_component(const vector < string > & args, Circuit & circuit) {
        if (args.size() != 1) throw runtime_error("Usage: delete <comp_name>");
        circuit.delete_component(args[0]);
    }
    void handle_rename_node(const vector < string > & args, Circuit & circuit) {
        if (args.size() != 3 || to_lower_util(args[0]) != "node") throw runtime_error("Invalid syntax");
        circuit.rename_node(args[1], args[2]);
    }
    void handle_add_component(const vector < string > & args, Circuit & circuit) {
        if (args.size() < 4) throw runtime_error("Insufficient arguments for add command.");
        string name = args[0];
        char type_char = toupper(name[0]);
        string n1 = args[1], n2 = args[2];
        vector < string > remaining_args(args.begin() + 3, args.end());

        switch (type_char) {
            case 'R':
            case 'C':
            case 'L':
                if (args.size() != 4) throw runtime_error("Incorrect argument count for " + string(1, type_char));
                if (type_char == 'R') circuit.add_component(make_unique < Resistor > (name, n1, n2, args[3]));
                else if (type_char == 'C') circuit.add_component(make_unique < Capacitor > (name, n1, n2, args[3]));
                else if (type_char == 'L') circuit.add_component(make_unique < Inductor > (name, n1, n2, args[3]));
                break;
            case 'I':
                if (args.size() != 4) throw runtime_error("Incorrect argument count for I source");
                circuit.add_component(make_unique < CurrentSource > (name, n1, n2, args[3]));
                break;
            case 'V':
                circuit.add_component(make_unique < VoltageSource > (name, n1, n2, remaining_args));
                break;
            case 'D':
                if (args.size() != 4) throw runtime_error("Incorrect argument count for Diode (D). Usage: D<name> n1 n2 ideal");
                circuit.add_component(make_unique < Diode > (name, n1, n2, args[3]));
                break;
            case 'E':
                if (args.size() != 6) throw runtime_error("Incorrect argument count for VCVS (E).");
                circuit.add_component(make_unique < VCVS > (name, n1, n2, args[3], args[4], args[5]));
                break;
            case 'G':
                if (args.size() != 6) throw runtime_error("Incorrect argument count for VCCS (G).");
                circuit.add_component(make_unique < VCCS > (name, n1, n2, args[3], args[4], args[5]));
                break;
            case 'H':
                if (args.size() != 5) throw runtime_error("Incorrect argument count for CCVS (H).");
                circuit.add_component(make_unique < CCVS > (name, n1, n2, args[3], args[4]));
                break;
            case 'F':
                if (args.size() != 5) throw runtime_error("Incorrect argument count for CCCS (F).");
                circuit.add_component(make_unique < CCCS > (name, n1, n2, args[3], args[4]));
                break;
            default:
                throw runtime_error("Element " + name + " not found in library");
        }
        cout << "Added component " << name << endl;
    }
    void handle_print(const vector < string > & args, Circuit & circuit) {
        if (args.size() < 2) throw runtime_error("Syntax error. Expected: .print <analysis_type> ...");
        string analysis_type = to_lower_util(args[0]);
        vector < string > vars_to_print;

        auto validate_vars = [&](const vector<string>& vars) {
            set<string> all_nodes;
            for (const auto& comp : circuit.components) {
                all_nodes.insert(comp->node1_name);
                all_nodes.insert(comp->node2_name);
            }
            if (circuit.ground_node_exists) {
                all_nodes.insert(circuit.ground_node_explicit_name);
                if (circuit.is_ground("0")) all_nodes.insert("0");
                if (circuit.is_ground("GND")) all_nodes.insert("GND");
            }

            for (const auto& var : vars) {
                string lower_var = to_lower_util(var);
                if (lower_var.rfind("v(", 0) == 0 && lower_var.back() == ')') {
                    string node_name = var.substr(2, var.length() - 3);
                    if (all_nodes.find(node_name) == all_nodes.end()) {
                        throw runtime_error("Node <" + node_name + "> not found in circuit");
                    }
                } else if (lower_var.rfind("i(", 0) == 0 && lower_var.back() == ')') {
                    string comp_name = var.substr(2, var.length() - 3);
                    auto it = find_if(circuit.components.begin(), circuit.components.end(), [&](const auto& comp){
                        return comp->name == comp_name;
                    });
                    if (it == circuit.components.end()) {
                        throw runtime_error("Component <" + comp_name + "> not found in circuit");
                    }
                } else {
                    throw runtime_error("Invalid print variable format: " + var + ". Use V(node) or I(component).");
                }
            }
        };

        if (analysis_type == "tran") {
            if (args.size() < 4) throw runtime_error("Syntax error. Expected: .print tran <Tstep> <Tstop> <var1> ...");
            double t_step = parse_value_with_metric_prefix_util(args[1]), t_stop = parse_value_with_metric_prefix_util(args[2]);
            vars_to_print.assign(args.begin() + 3, args.end());

            validate_vars(vars_to_print);

            circuit.perform_transient_analysis(t_step, t_stop);
            if (!circuit.tran_solved) return;
            cout << left << setw(15) << "time (s)";
            for (const auto &
                        var: vars_to_print) cout << setw(15) <<
                                                      var;
            cout << endl;
            for (const auto & result_point: circuit.tran_results) {
                cout << fixed << setprecision(6) << setw(15) << result_point.at("time");
                for (const auto &
                            var: vars_to_print) {
                    try {
                        ostringstream val_ss;
                        val_ss << fixed << setprecision(6) << result_point.at(var);
                        string val_str = val_ss.str();
                        char type = to_lower_util(var)[0];
                        if (type == 'v') val_str += " V";
                        else if (type == 'i') val_str += " A";
                        cout << setw(15) << val_str;
                    } catch (const out_of_range & ) {
                        cout << setw(15) << "N/A";
                    }
                }
                cout << endl;
            }
        } else if (analysis_type == "dc") {
            if (args.size() < 6) throw runtime_error("Syntax error. Expected: .print dc <SrcName> <Start> <Stop> <Incr> <var1> ...");
            string src_name = args[1];
            double start_val = parse_value_with_metric_prefix_util(args[2]), end_val = parse_value_with_metric_prefix_util(args[3]), increment = parse_value_with_metric_prefix_util(args[4]);
            vars_to_print.assign(args.begin() + 5, args.end());

            validate_vars(vars_to_print);

            circuit.perform_dc_sweep_analysis(src_name, start_val, end_val, increment);
            if (!circuit.dc_sweep_solved) return;

            string header_src_name = circuit.dc_sweep_source_name;
            char src_name_char = to_lower_util(header_src_name)[0];
            if (src_name_char == 'v') header_src_name += " (V)";
            else if (src_name_char == 'i') header_src_name += " (A)";
            cout << left << setw(15) << header_src_name;

            for (const auto &
                        var: vars_to_print) cout << setw(15) <<
                                                      var;
            cout << endl;
            for (const auto & result_point: circuit.dc_sweep_results) {
                cout << fixed << setprecision(6) << setw(15) << result_point.at(circuit.dc_sweep_source_name);
                for (const auto &
                            var: vars_to_print) {
                    if (result_point.count(var)) {
                        ostringstream val_ss;
                        val_ss << fixed << setprecision(6) << result_point.at(var);
                        string val_str = val_ss.str();
                        char type = to_lower_util(var)[0];
                        if (type == 'v') val_str += " V";
                        else if (type == 'i') val_str += " A";
                        cout << setw(15) << val_str;
                    } else {
                        cout << setw(15) << "N/A";
                    }
                }
                cout << endl;
            }
        } else throw runtime_error("Unsupported analysis type '" + analysis_type + "'.");
    }
    void handle_list_nodes(const Circuit & circuit) {
        set < string > all_nodes;
        for (const auto & comp: circuit.components) {
            all_nodes.insert(comp -> node1_name);
            all_nodes.insert(comp -> node2_name);
        }
        cout << "Available nodes:" << endl;
        if (all_nodes.empty()) {
            cout << "  (No components in circuit)" << endl;
            return;
        }
        for (const auto & node_name: all_nodes) cout << "  " << node_name << (circuit.is_ground(node_name) ? " (Ground)" : "") << endl;
    }
    void handle_list_components(const vector < string > & args,
                                const Circuit & circuit) {
        string filter_type_str;
        if (!args.empty()) filter_type_str = to_lower_util(args[0]);
        cout << "Components in circuit" << (filter_type_str.empty() ? "" : " (filtered by type: " + filter_type_str + ")") << ":" << endl;
        if (circuit.components.empty()) {
            cout << "  (No components)" << endl;
            return;
        }
        bool found_any = false;
        for (const auto & comp: circuit.components) {
            char comp_type_char = to_lower_util(comp -> name)[0];
            if (filter_type_str.empty() || filter_type_str.find(comp_type_char) != string::npos) {
                found_any = true;
                cout << "  - " << comp -> to_netlist_string() << endl;
            }
        }
        if (!found_any) cout << "  (No components of specified type '" << filter_type_str << "' found)" << endl;
    }
};

void print_banner() {
    cout << "---------------------------------------------------" << endl;
    cout << "                     SUTSpice" << endl;
    cout << "           Erfan Rahmati - Ilia Ghaderi" << endl;
    cout << "    OOP 2025 - Sharif University of Technology" << endl;
    cout << "---------------------------------------------------" << endl;
    cout << "Circuit Commands: " << endl;
    cout << "  add <comp> ..." << endl;
    cout << "  delete <comp_name>" << endl;
    cout << "  .rename node <old> <new>" << endl;
    cout << "  .nodes" << endl;
    cout << "  .list [type]" << endl;
    cout << "\nFile Commands:" << endl;
    cout << "  .load <filepath>" << endl;
    cout << "  .show existing schematics" << endl;
    cout << "  .save [filename]" << endl;
    cout << "  .clear" << endl;
    cout << "\nAnalysis Commands:" << endl;
    cout << "  .print tran <Tstep> <Tstop> <...>" << endl;
    cout << "  .print dc <Src> <Start> <End> <Inc> ..." << endl;
    cout << "---------------------------------------------------" << endl;
}

int main() {
    Circuit circuit_main;
    Parser parser_main;
    string line_input;
    print_banner();
    while (true) {
        if (!getline(cin, line_input)) break;
        line_input = trim_string_util(line_input);
        if (line_input.empty()) continue;

        Command cmd_main = parser_main.parse_line(line_input);
        string cmd_type_lower = to_lower_util(cmd_main.type);

        if (cmd_type_lower == "exit" || cmd_type_lower == "quit") {
            if (circuit_main.is_dirty) {
                cout << "You have unsaved changes. Save before exiting? (yes/no): ";
                string response;
                getline(cin, response);
                if (to_lower_util(trim_string_util(response)) == "yes") {
                    if (!circuit_main.current_filepath.empty()) {
                        try {
                            circuit_main.save_to_file(circuit_main.current_filepath);
                        } catch (const exception& e) {
                            cerr << "Error saving file: " << e.what() << endl;
                        }
                    } else {
                        parser_main.execute_command({".save", {}}, circuit_main);
                    }
                }
            }
            cout << "Exiting simulator." << endl;
            break;
        }
        parser_main.execute_command(cmd_main, circuit_main);
    }
    return 0;
}
