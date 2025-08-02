#include <SDL2/SDL.h>
#include <SDL2/SDL_image.h>
#include <SDL2/SDL_ttf.h>
#ifdef _WIN32
#include <SDL2/SDL2_gfx.h>
#else
#include <SDL2/SDL2_gfxPrimitives.h>
#endif

#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <functional>
#include <map>
#include <set>
#include <memory>
#include <stdexcept>
#include <algorithm>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <numeric>

// Windows specific for directory listing
#ifdef _WIN32
#include <windows.h>
#include <direct.h>
#endif

using namespace std;

// --- Constants ---
const int SCREEN_WIDTH = 1280;
const int SCREEN_HEIGHT = 720;
const int TOP_BAR_HEIGHT = 40;
const int GRID_SPACING = 20;
const int COMPONENT_DEFAULT_LENGTH = GRID_SPACING * 4;

const string BASE_PATH = "C:/Users/Erfan/Dev/Cpp/sutSpice_phase2/";
const string ASSET_PATH = BASE_PATH + "assets/";
const string SCHEMATICS_PATH = BASE_PATH + "schematics/";

// --- Utility Functions (Shared) ---
inline string to_lower_util(string s) {
    transform(s.begin(), s.end(), s.begin(),
              [](unsigned char c) { return tolower(c); });
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
        base_val = stod(val_str, &suffix_start_pos);
    }
    catch (const invalid_argument&) {
        throw invalid_argument("Invalid numeric value in string: \"" + val_str_orig + "\"");
    }
    catch (const out_of_range&) {
        throw out_of_range("Numeric value out of range in string: \"" + val_str_orig + "\"");
    }
    if (suffix_start_pos >= val_str.length()) {
        return base_val;
    }
    string suffix_str = to_lower_util(val_str.substr(suffix_start_pos));
    double multiplier = 1.0;
    if (suffix_str.rfind("meg", 0) == 0) {
        multiplier = 1e6;
    }
    else if (!suffix_str.empty()) {
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

vector<double> gaussian_elimination_matrix(vector<vector<double>> A, vector<double> b) {
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
    vector<double> x(n);
    for (int i = n - 1; i >= 0; --i) {
        if (abs(A[i][i]) < 1e-12) {
            if (abs(b[i]) > 1e-12) {
                throw runtime_error("System is inconsistent (no solution).");
            }
            x[i] = 0;
        }
        else {
            x[i] = b[i];
            for (int j = i + 1; j < n; ++j) x[i] -= A[i][j] * x[j];
            x[i] /= A[i][i];
        }
    }
    return x;
}

double dist_to_segment_sq(SDL_Point p, SDL_Point v, SDL_Point w) {
    double l2 = pow(v.x - w.x, 2) + pow(v.y - w.y, 2);
    if (l2 == 0.0) return pow(p.x - v.x, 2) + pow(p.y - v.y, 2);
    double t = max(0.0, min(1.0, ((p.x - v.x) * (w.x - v.x) + (p.y - v.y) * (w.y - v.y)) / l2));
    return pow(p.x - (v.x + t * (w.x - v.x)), 2) + pow(p.y - (v.y + t * (w.y - v.y)), 2);
}


// --- SPICE Simulation Engine (from first program, wrapped in a namespace) ---
namespace SpiceEngine {

    class Circuit;

    enum class ComponentType {
        Resistor, Capacitor, Inductor, VoltageSource, CurrentSource,
        VCVS, VCCS, CCVS, CCCS, Diode
    };
    using ResultPoint = map<string, double>;
    const string GROUND_NODE_NAME_CONST = "0";

    class Component {
    public:
        string name;
        string node1_name, node2_name;
        double value;
        Component(string name_val, string n1, string n2, double val = 0.0) : name(move(name_val)),
                                                                             node1_name(move(n1)), node2_name(move(n2)), value(val) {}
        Component(string name_val, string n1, string n2, const string& val_str) : name(move(name_val)),
                                                                                  node1_name(move(n1)), node2_name(move(n2)) {
            this->value = parse_value_with_metric_prefix_util(val_str);
        }
        virtual ~Component() = default;
        virtual ComponentType get_type() const = 0;
        virtual string to_netlist_string() const = 0;
        virtual void stamp(Circuit&, vector<vector<double>>&, vector<vector<double>>&, vector<vector<double>>&, vector<vector<double>>&, vector<double>&, vector<double>&, map<string, int>&, double,
                           const vector<double>&) = 0;
        virtual void update_time_dependant_value(double) {}
    };

    class Resistor : public Component {
    public:
        Resistor(const string& r_name, const string& n1, const string& n2, const string& val_str) : Component(r_name, n1, n2, val_str) {
            if (value <= 0) throw runtime_error("Resistor " + name + " must have positive resistance.");
        }
        ComponentType get_type() const override { return ComponentType::Resistor; }
        string to_netlist_string() const override {
            ostringstream oss;
            oss << name << " " << node1_name << " " << node2_name << " " << fixed << setprecision(12) << value;
            return oss.str();
        }
        void stamp(Circuit& circuit, vector<vector<double>>& G, vector<vector<double>>&, vector<vector<double>>&, vector<vector<double>>&, vector<double>&, vector<double>&, map<string, int>&, double, const vector<double>&) override;
    };

    class VoltageSource : public Component {
    public:
        enum SourceType { DC, SINUSOIDAL, PULSE };
        SourceType sourceType;
        double dc_offset, amplitude, frequency;
        double v1, v2, td, tr, tf, pw, per;
        vector<string> raw_params;
        VoltageSource(const string& v_name, const string& n1, const string& n2, const vector<string>& params);
        ComponentType get_type() const override { return ComponentType::VoltageSource; }
        string to_netlist_string() const override;
        void stamp(Circuit& circuit, vector<vector<double>>&, vector<vector<double>>& B, vector<vector<double>>& C, vector<vector<double>>&, vector<double>&, vector<double>& E, map<string, int>& m_map, double, const vector<double>&) override;
        void update_time_dependant_value(double time) override;
    };

    class Capacitor : public Component {
    public:
        Capacitor(const string& c_name, const string& n1, const string& n2, const string& val_str) : Component(c_name, n1, n2, val_str) {
            if (value <= 0) throw runtime_error("Capacitance must be positive.");
        }
        ComponentType get_type() const override { return ComponentType::Capacitor; }
        string to_netlist_string() const override {
            ostringstream oss;
            oss << name << " " << node1_name << " " << node2_name << " " << fixed << setprecision(12) << value;
            return oss.str();
        }
        void stamp(Circuit& circuit, vector<vector<double>>& G, vector<vector<double>>&, vector<vector<double>>&, vector<vector<double>>&, vector<double>& J, vector<double>&, map<string, int>&, double h, const vector<double>& prev_sol) override;
    };

    class Inductor : public Component {
    public:
        Inductor(const string& l_name, const string& n1, const string& n2, const string& val_str) : Component(l_name, n1, n2, val_str) {
            if (value <= 0) throw runtime_error("Inductance must be positive.");
        }
        ComponentType get_type() const override { return ComponentType::Inductor; }
        string to_netlist_string() const override {
            ostringstream oss;
            oss << name << " " << node1_name << " " << node2_name << " " << fixed << setprecision(12) << value;
            return oss.str();
        }
        void stamp(Circuit& circuit, vector<vector<double>>&, vector<vector<double>>& B, vector<vector<double>>& C, vector<vector<double>>& D, vector<double>&, vector<double>& E, map<string, int>& m_map, double h, const vector<double>& prev_sol) override;
    };

    class CurrentSource : public Component {
    public:
        CurrentSource(const string& i_name, const string& n1, const string& n2, const string& val_str) : Component(i_name, n1, n2, val_str) {}
        ComponentType get_type() const override { return ComponentType::CurrentSource; }
        string to_netlist_string() const override {
            ostringstream oss;
            oss << name << " " << node1_name << " " << node2_name << " " << fixed << setprecision(12) << value;
            return oss.str();
        }
        void stamp(Circuit& circuit, vector<vector<double>>&, vector<vector<double>>&, vector<vector<double>>&, vector<vector<double>>&, vector<double>& J, vector<double>&, map<string, int>&, double, const vector<double>&) override;
    };

    class Diode : public Component {
    public:
        bool is_on = false;
        const double Ron = 1e-3;
        const double Roff = 1e9;
        Diode(const string& name, const string& n1, const string& n2, const string& model) : Component(name, n1, n2, 0.0) {
            if (to_lower_util(model) != "ideal" && to_lower_util(model) != "1n4148") { // Allow both models
                throw runtime_error("Only 'ideal' or '1N4148' diode model is supported. Got: " + model);
            }
        }
        ComponentType get_type() const override { return ComponentType::Diode; }
        string to_netlist_string() const override { return name + " " + node1_name + " " + node2_name + " ideal"; }
        void stamp(Circuit& circuit, vector<vector<double>>& G, vector<vector<double>>&, vector<vector<double>>&, vector<vector<double>>&, vector<double>&, vector<double>&, map<string, int>&, double, const vector<double>&) override;
    };

    class Circuit {
    public:
        vector<unique_ptr<Component>> components;
        map<string, int> node_to_idx;
        vector<string> idx_to_node_name;
        string ground_node_explicit_name = GROUND_NODE_NAME_CONST;
        bool ground_node_exists = false;
        vector<VoltageSource*> voltage_source_list;
        vector<Inductor*> inductor_list;
        vector<Diode*> diode_list;
        vector<Component*> vcvs_list;
        vector<Component*> ccvs_list;
        bool tran_solved = false;
        vector<ResultPoint> tran_results;
        bool dc_sweep_solved = false;
        vector<ResultPoint> dc_sweep_results;
        string dc_sweep_source_name;

        void set_ground_node(const string& node_name) {
            ground_node_explicit_name = node_name;
            ground_node_exists = true;
        }
        bool is_ground(const string& node_name) const {
            if (!ground_node_exists) return false;
            if (node_name == ground_node_explicit_name) return true;
            if (ground_node_explicit_name == "0" && node_name == "GND") return true;
            if (ground_node_explicit_name == "GND" && node_name == "0") return true;
            return false;
        }
        void add_component(unique_ptr<Component> comp) {
            for (const auto& existing_comp : components) {
                if (existing_comp->name == comp->name) {
                    throw runtime_error("Error: Component with name '" + comp->name + "' already exists.");
                }
            }
            components.push_back(move(comp));
        }
        int prepare_for_analysis();
        int get_node_matrix_index(const string& node_name) const;
        double get_voltage_at(const string& node_name, const ResultPoint& results) const;
        void build_mna_matrix(vector<vector<double>>& A, vector<double>& z, double h, const vector<double>& prev_sol);
        void calculate_and_store_passive_currents(ResultPoint& result_point, const ResultPoint& prev_result_point, double h);
        void perform_transient_analysis(double t_step, double t_stop);
        void perform_dc_sweep_analysis(const string& src, double start, double end, double inc);
    };

    void Resistor::stamp(Circuit& circuit, vector<vector<double>>& G, vector<vector<double>>&, vector<vector<double>>&, vector<vector<double>>&, vector<double>&, vector<double>&, map<string, int>&, double, const vector<double>&) {
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

    VoltageSource::VoltageSource(const string& v_name, const string& n1, const string& n2, const vector<string>& params)
            : Component(v_name, n1, n2, 0.0), sourceType(DC), dc_offset(0.0), amplitude(0.0), frequency(0.0), v1(0), v2(0), td(0), tr(0), tf(0), pw(0), per(0), raw_params(params) {
        if (params.empty()) throw runtime_error("No value or parameters provided for voltage source " + name);

        string combined_params;
        for (size_t i = 0; i < params.size(); ++i) {
            combined_params += params[i] + (i < params.size() - 1 ? " " : "");
        }

        string first_param_upper = to_lower_util(params[0]);
        if (first_param_upper.rfind("sin", 0) == 0) {
            size_t start_paren = combined_params.find('('), end_paren = combined_params.rfind(')');
            if (start_paren == string::npos || end_paren == string::npos) throw runtime_error("Mismatched parentheses in SIN() for " + name);
            stringstream ss(combined_params.substr(start_paren + 1, end_paren - start_paren - 1));
            string offset_str, amp_str, freq_str;
            ss >> offset_str >> amp_str >> freq_str;
            if (ss.fail() || !ss.eof()) throw runtime_error("Invalid parameters inside SIN() for " + name);
            sourceType = SINUSOIDAL;
            dc_offset = parse_value_with_metric_prefix_util(offset_str);
            amplitude = parse_value_with_metric_prefix_util(amp_str);
            frequency = parse_value_with_metric_prefix_util(freq_str);
            this->value = dc_offset;
        }
        else if (first_param_upper.rfind("pulse", 0) == 0) {
            size_t start_paren = combined_params.find('('), end_paren = combined_params.rfind(')');
            if (start_paren == string::npos || end_paren == string::npos) throw runtime_error("Mismatched parentheses in PULSE() for " + name);
            stringstream ss(combined_params.substr(start_paren + 1, end_paren - start_paren - 1));
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
            this->value = v1;
        }
        else if (params.size() == 1) {
            sourceType = DC;
            this->value = parse_value_with_metric_prefix_util(params[0]);
            this->dc_offset = this->value;
        }
        else throw runtime_error("Invalid parameters for voltage source " + name);
    }

    string VoltageSource::to_netlist_string() const {
        ostringstream oss;
        oss << name << " " << node1_name << " " << node2_name << " ";
        for (size_t i = 0; i < raw_params.size(); ++i) {
            oss << raw_params[i] << (i == raw_params.size() - 1 ? "" : " ");
        }
        return oss.str();
    }

    void VoltageSource::stamp(Circuit& circuit, vector<vector<double>>&, vector<vector<double>>& B, vector<vector<double>>& C, vector<vector<double>>&, vector<double>&, vector<double>& E, map<string, int>& m_map, double, const vector<double>&) {
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
        E[m_idx] += this->value;
    }

    void VoltageSource::update_time_dependant_value(double time) {
        if (sourceType == SINUSOIDAL) {
            value = (frequency > 0) ? (dc_offset + amplitude * sin(2 * 3.14159265358979323846 * frequency * time)) : (dc_offset + amplitude);
        }
        else if (sourceType == PULSE) {
            if (time < td) value = v1;
            else {
                double t_rel = per > 0 ? fmod(time - td, per) : time - td;
                if (tr > 0 && t_rel < tr) value = v1 + (v2 - v1) * (t_rel / tr);
                else if (t_rel < tr + pw) value = v2;
                else if (tf > 0 && t_rel < tr + pw + tf) value = v2 + (v1 - v2) * ((t_rel - tr - pw) / tf);
                else value = v1;
            }
        }
    }

    void CurrentSource::stamp(Circuit& circuit, vector<vector<double>>&, vector<vector<double>>&, vector<vector<double>>&, vector<vector<double>>&, vector<double>& J, vector<double>&, map<string, int>&, double, const vector<double>&) {
        int idx1 = circuit.get_node_matrix_index(node1_name);
        int idx2 = circuit.get_node_matrix_index(node2_name);
        if (idx1 >= 0) J[idx1] -= value;
        if (idx2 >= 0) J[idx2] += value;
    }

    void Capacitor::stamp(Circuit& circuit, vector<vector<double>>& G, vector<vector<double>>&, vector<vector<double>>&, vector<vector<double>>&, vector<double>& J, vector<double>&, map<string, int>&, double h, const vector<double>& prev_sol) {
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

    void Inductor::stamp(Circuit& circuit, vector<vector<double>>&, vector<vector<double>>& B, vector<vector<double>>& C, vector<vector<double>>& D, vector<double>&, vector<double>& E, map<string, int>& m_map, double h, const vector<double>& prev_sol) {
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

    void Diode::stamp(Circuit& circuit, vector<vector<double>>& G, vector<vector<double>>&, vector<vector<double>>&, vector<vector<double>>&, vector<double>&, vector<double>&, map<string, int>&, double, const vector<double>&) {
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

    int Circuit::prepare_for_analysis() {
        node_to_idx.clear();
        idx_to_node_name.clear();
        voltage_source_list.clear();
        inductor_list.clear();
        vcvs_list.clear();
        ccvs_list.clear();
        diode_list.clear();
        set<string> unique_node_names;
        for (const auto& comp : components) {
            unique_node_names.insert(comp->node1_name);
            unique_node_names.insert(comp->node2_name);
            if (auto vs = dynamic_cast<VoltageSource*>(comp.get())) voltage_source_list.push_back(vs);
            else if (auto ind = dynamic_cast<Inductor*>(comp.get())) inductor_list.push_back(ind);
            else if (auto d = dynamic_cast<Diode*>(comp.get())) diode_list.push_back(d);
        }
        if (!ground_node_exists && (unique_node_names.count("0") || unique_node_names.count("GND"))) {
            set_ground_node(unique_node_names.count("0") ? "0" : "GND");
        }
        int idx = 0;
        for (const auto& name : unique_node_names)
            if (!is_ground(name)) {
                node_to_idx[name] = idx++;
                idx_to_node_name.push_back(name);
            }
        return idx;
    }

    int Circuit::get_node_matrix_index(const string& node_name) const {
        if (is_ground(node_name)) return -1;
        auto it = node_to_idx.find(node_name);
        if (it != node_to_idx.end()) {
            return it->second;
        }
        return -2;
    }

    double Circuit::get_voltage_at(const string& node_name, const ResultPoint& results) const {
        if (is_ground(node_name)) return 0.0;
        string v_name = "V(" + node_name + ")";
        auto it = results.find(v_name);
        if (it != results.end()) return it->second;
        return 0.0;
    }

    void Circuit::build_mna_matrix(vector<vector<double>>& A, vector<double>& z, double h, const vector<double>& prev_sol) {
        if (!components.empty() && !ground_node_exists) throw runtime_error("No ground node defined.");
        int N = idx_to_node_name.size();
        int M = voltage_source_list.size() + inductor_list.size() + vcvs_list.size() + ccvs_list.size();
        int system_size = N + M;
        if (system_size == 0 && diode_list.empty()) {
            A.clear(); z.clear(); return;
        }
        A.assign(system_size, vector<double>(system_size, 0.0));
        z.assign(system_size, 0.0);
        vector<vector<double>> G(N, vector<double>(N, 0.0)), B(N, vector<double>(M, 0.0)), C(M, vector<double>(N, 0.0)), D(M, vector<double>(M, 0.0));
        vector<double> J(N, 0.0), E(M, 0.0);
        map<string, int> m_map;
        int m_counter = 0;
        for (const auto& vs : voltage_source_list) m_map[vs->name] = m_counter++;
        for (const auto& l : inductor_list) m_map[l->name] = m_counter++;
        for (const auto& comp : components) comp->stamp(*this, G, B, C, D, J, E, m_map, h, prev_sol);
        for (int i = 0; i < N; ++i) for (int j = 0; j < N; ++j) A[i][j] = G[i][j];
        for (int i = 0; i < N; ++i) for (int j = 0; j < M; ++j) A[i][N + j] = B[i][j];
        for (int i = 0; i < M; ++i) for (int j = 0; j < N; ++j) A[N + i][j] = C[i][j];
        for (int i = 0; i < M; ++i) for (int j = 0; j < M; ++j) A[N + i][N + j] = D[i][j];
        for (int i = 0; i < N; ++i) z[i] = J[i];
        for (int i = 0; i < M; ++i) z[N + i] = E[i];
    }

    void Circuit::calculate_and_store_passive_currents(ResultPoint& result_point, const ResultPoint& prev_result_point, double h) {
        for (const auto& comp : components) {
            double v1 = get_voltage_at(comp->node1_name, result_point);
            double v2 = get_voltage_at(comp->node2_name, result_point);
            if (auto r = dynamic_cast<Resistor*>(comp.get())) {
                result_point["I(" + r->name + ")"] = (v1 - v2) / r->value;
            }
            else if (auto c = dynamic_cast<Capacitor*>(comp.get())) {
                if (h > 0 && !prev_result_point.empty()) {
                    double v1_prev = get_voltage_at(c->node1_name, prev_result_point);
                    double v2_prev = get_voltage_at(c->node2_name, prev_result_point);
                    result_point["I(" + c->name + ")"] = c->value * ((v1 - v2) - (v1_prev - v2_prev)) / h;
                }
                else {
                    result_point["I(" + c->name + ")"] = 0.0;
                }
            }
            else if (auto d = dynamic_cast<Diode*>(comp.get())) {
                result_point["I(" + d->name + ")"] = (v1 - v2) / (d->is_on ? d->Ron : d->Roff);
            }
        }
    }

    void Circuit::perform_transient_analysis(double t_step, double t_stop) {
        tran_solved = false;
        tran_results.clear();
        if (t_step <= 0 || t_stop <= 0 || t_step > t_stop) throw runtime_error("Invalid transient parameters.");
        int N = prepare_for_analysis();
        int M = voltage_source_list.size() + inductor_list.size() + vcvs_list.size() + ccvs_list.size();
        if (N + M == 0 && diode_list.empty()) {
            tran_solved = true;
            return;
        }
        for (auto* d : diode_list) { d->is_on = false; }
        vector<double> prev_solution(N + M, 0.0);
        ResultPoint prev_result_point;
        for (double t = 0; t <= t_stop + (t_step / 2.0); t += t_step) {
            for (auto& comp : components) { comp->update_time_dependant_value(t); }
            const int MAX_DIODE_ITERATIONS = 100;
            bool diodes_converged = false;
            vector<double> current_solution;
            for (int iter = 0; iter < MAX_DIODE_ITERATIONS; ++iter) {
                vector<vector<double>> A;
                vector<double> z;
                build_mna_matrix(A, z, t_step, prev_solution);
                current_solution = gaussian_elimination_matrix(A, z);
                if (diode_list.empty()) { diodes_converged = true; break; }
                bool state_changed = false;
                for (auto* d : diode_list) {
                    bool old_state = d->is_on;
                    double v1 = (get_node_matrix_index(d->node1_name) >= 0) ? current_solution[get_node_matrix_index(d->node1_name)] : 0.0;
                    double v2 = (get_node_matrix_index(d->node2_name) >= 0) ? current_solution[get_node_matrix_index(d->node2_name)] : 0.0;
                    d->is_on = (v1 > v2);
                    if (d->is_on != old_state) { state_changed = true; }
                }
                if (!state_changed) { diodes_converged = true; break; }
            }
            if (!diodes_converged) { throw runtime_error("Diode model failed to converge at time " + to_string(t)); }
            ResultPoint result_at_t;
            result_at_t["time"] = t;
            for (int i = 0; i < N; ++i) result_at_t["V(" + idx_to_node_name[i] + ")"] = current_solution[i];
            if (ground_node_exists) result_at_t["V(" + ground_node_explicit_name + ")"] = 0.0;
            map<string, int> m_map;
            int m_counter = 0;
            for (const auto& vs : voltage_source_list) m_map[vs->name] = m_counter++;
            for (const auto& l : inductor_list) m_map[l->name] = m_counter++;
            for (const auto& p : m_map) result_at_t["I(" + p.first + ")"] = current_solution[N + p.second];
            calculate_and_store_passive_currents(result_at_t, prev_result_point, t_step);
            tran_results.push_back(result_at_t);
            prev_solution = current_solution;
            prev_result_point = result_at_t;
        }
        tran_solved = true;
    }

    void Circuit::perform_dc_sweep_analysis(const string& src, double start, double end, double inc) {
        dc_sweep_solved = false;
        dc_sweep_results.clear();
        if (inc == 0 || (end > start && inc < 0) || (end < start && inc > 0)) throw runtime_error("Invalid sweep parameters.");
        Component* sweep_comp = nullptr;
        for (auto& c : components)
            if (c->name == src) {
                sweep_comp = c.get();
                break;
            }
        if (!sweep_comp || (sweep_comp->get_type() != ComponentType::VoltageSource && sweep_comp->get_type() != ComponentType::CurrentSource)) throw runtime_error("Sweep source not found or invalid.");
        double orig_val = sweep_comp->value;
        dc_sweep_source_name = src;
        int N = prepare_for_analysis();
        int M = voltage_source_list.size() + inductor_list.size() + vcvs_list.size() + ccvs_list.size();
        if (N + M == 0 && diode_list.empty()) {
            dc_sweep_solved = true;
            return;
        }
        for (auto* d : diode_list) { d->is_on = false; }
        for (double val = start; (inc > 0) ? (val <= end + abs(inc) / 2.0) : (val >= end - abs(inc) / 2.0); val += inc) {
            sweep_comp->value = val;
            const int MAX_DIODE_ITERATIONS = 100;
            bool diodes_converged = false;
            vector<double> solution;
            for (int iter = 0; iter < MAX_DIODE_ITERATIONS; ++iter) {
                vector<vector<double>> A;
                vector<double> z;
                build_mna_matrix(A, z, 0.0, {});
                solution = gaussian_elimination_matrix(A, z);
                if (diode_list.empty()) { diodes_converged = true; break; }
                bool state_changed = false;
                for (auto* d : diode_list) {
                    bool old_state = d->is_on;
                    double v1 = (get_node_matrix_index(d->node1_name) >= 0 && get_node_matrix_index(d->node1_name) < solution.size()) ? solution[get_node_matrix_index(d->node1_name)] : 0.0;
                    double v2 = (get_node_matrix_index(d->node2_name) >= 0 && get_node_matrix_index(d->node2_name) < solution.size()) ? solution[get_node_matrix_index(d->node2_name)] : 0.0;
                    d->is_on = (v1 > v2);
                    if (d->is_on != old_state) { state_changed = true; }
                }
                if (!state_changed) { diodes_converged = true; break; }
            }
            if (!diodes_converged) { throw runtime_error("Diode model failed to converge at sweep value " + to_string(val)); }
            ResultPoint result_at_val;
            result_at_val[src] = val;
            for (int i = 0; i < N; ++i) result_at_val["V(" + idx_to_node_name[i] + ")"] = solution[i];
            if (ground_node_exists) result_at_val["V(" + ground_node_explicit_name + ")"] = 0.0;
            map<string, int> m_map;
            int m_counter = 0;
            for (const auto& vs : voltage_source_list) m_map[vs->name] = m_counter++;
            for (const auto& l : inductor_list) m_map[l->name] = m_counter++;
            for (const auto& p : m_map) result_at_val["I(" + p.first + ")"] = solution[N + p.second];
            calculate_and_store_passive_currents(result_at_val, {}, 0.0);
            dc_sweep_results.push_back(result_at_val);
        }
        sweep_comp->value = orig_val;
        dc_sweep_solved = true;
    }

} // end namespace SpiceEngine

// --- GUI Data Structures ---
bool isFileMenuOpen = false;
bool isComponentsMenuOpen = false;
bool isSimulateMenuOpen = false;
enum class InteractionMode {
    NONE, WIRING, PLACE_COMPONENT, DRAGGING_COMPONENT, DELETE_ITEM,
    EDITING_COMPONENT_VALUE, PLACE_LABEL, EDITING_LABEL_TEXT, PLACE_GND_LABEL, DIALOG_ACTIVE
};
InteractionMode currentInteractionMode = InteractionMode::NONE;
bool isDrawingWire = false;
SDL_Point firstInteractionPoint = { -1, -1 };

enum class ComponentType {
    NONE, RESISTOR, CAPACITOR, INDUCTOR, DIODE,
    DC_VOLTAGE_SOURCE, DC_CURRENT_SOURCE, AC_VOLTAGE_SOURCE
};
ComponentType selectedComponentType = ComponentType::NONE;
int placementRotation = 0;

struct Wire { SDL_Point start, end; };
struct Component {
    ComponentType type;
    SDL_Point node1, node2;
    string id, value;
    int rotation_angle = 0;
    SDL_Rect value_rect, id_rect;
};
struct NodeLabel {
    string text;
    SDL_Point position;
    SDL_Rect text_rect;
    int rotation_angle = 0;
};
vector<Wire> wires;
vector<Component> components;
vector<NodeLabel> labels;
Component* editingComponent = nullptr;
NodeLabel* editingLabel = nullptr;
string textInputBuffer = "";
string currentSchematicFileName = "";
bool schematicModified = false;
map<ComponentType, int> componentNameCounters;

enum class DialogType { NONE, SAVE_AS, OPEN, TRANSIENT_ANALYSIS, DC_SWEEP_ANALYSIS };
DialogType activeDialogType = DialogType::NONE;
struct DialogField {
    string label, buffer;
    SDL_Rect input_rect;
};
vector<DialogField> dialogFields;
string dialogTitle;
int activeDialogFieldIndex = -1;

class Button {
public:
    Button(string text, int x, int y, int w, int h, function<void()> onClick)
            : m_text(text), m_position({ x, y, w, h }), m_onClick(onClick), m_is_hovered(false) {}
    SDL_Rect m_position;
    function<void()> m_onClick;
    bool m_is_hovered;
    void handle_event(SDL_Event* e) {
        if (e->type == SDL_MOUSEBUTTONDOWN) {
            int x, y;
            SDL_GetMouseState(&x, &y);
            if (SDL_PointInRect(new const SDL_Point{ x, y }, &m_position)) m_onClick();
        }
        else if (e->type == SDL_MOUSEMOTION) {
            int x, y;
            SDL_GetMouseState(&x, &y);
            m_is_hovered = SDL_PointInRect(new const SDL_Point{ x, y }, &m_position);
        }
    }
    void render(SDL_Renderer* renderer, TTF_Font* font) {
        SDL_SetRenderDrawColor(renderer, m_is_hovered ? 0x66 : 0x44, m_is_hovered ? 0x66 : 0x44, m_is_hovered ? 0x66 : 0x44, 0xFF);
        SDL_RenderFillRect(renderer, &m_position);
        if (font && !m_text.empty()) {
            SDL_Color textColor = { 0xFF, 0xFF, 0xFF, 0xFF };
            SDL_Surface* textSurface = TTF_RenderText_Solid(font, m_text.c_str(), textColor);
            if (textSurface) {
                SDL_Texture* textTexture = SDL_CreateTextureFromSurface(renderer, textSurface);
                SDL_Rect renderQuad = { m_position.x + (m_position.w - textSurface->w) / 2, m_position.y + (m_position.h - textSurface->h) / 2, textSurface->w, textSurface->h };
                SDL_RenderCopy(renderer, textTexture, nullptr, &renderQuad);
                SDL_DestroyTexture(textTexture);
                SDL_FreeSurface(textSurface);
            }
        }
    }
private:
    string m_text;
};

struct ComponentMenuItem {
    string name;
    ComponentType type;
    SDL_Texture* iconTexture = nullptr;
    Button button;
};

// --- GUI Functions ---
SDL_Point snap_to_grid(int x, int y) {
    int snappedX = round((float)x / GRID_SPACING) * GRID_SPACING;
    int snappedY = round((float)y / GRID_SPACING) * GRID_SPACING;
    if (snappedY < TOP_BAR_HEIGHT) snappedY = TOP_BAR_HEIGHT;
    return { snappedX, snappedY };
}

void draw_grid(SDL_Renderer* renderer) {
    SDL_SetRenderDrawColor(renderer, 0x40, 0x40, 0x40, 0xFF);
    for (int x = 0; x < SCREEN_WIDTH; x += GRID_SPACING) {
        for (int y = TOP_BAR_HEIGHT; y < SCREEN_HEIGHT; y += GRID_SPACING) {
            SDL_RenderDrawPoint(renderer, x, y);
        }
    }
}

SDL_Texture* loadAndProcessTexture(const string& path, SDL_Renderer* renderer) {
    SDL_Surface* originalSurface = IMG_Load(path.c_str());
    if (!originalSurface) {
        cerr << "Failed to load image " << path << "! SDL_image Error: " << IMG_GetError() << endl;
        return nullptr;
    }
    SDL_LockSurface(originalSurface);
    auto pixels = static_cast<uint32_t*>(originalSurface->pixels);
    SDL_PixelFormat* format = originalSurface->format;
    uint32_t white = SDL_MapRGB(format, 255, 255, 255);
    for (int i = 0; i < originalSurface->w * originalSurface->h; ++i) {
        uint8_t r, g, b, a;
        SDL_GetRGBA(pixels[i], format, &r, &g, &b, &a);
        if (a > 50) pixels[i] = white;
    }
    SDL_UnlockSurface(originalSurface);
    SDL_Texture* finalTexture = SDL_CreateTextureFromSurface(renderer, originalSurface);
    SDL_FreeSurface(originalSurface);
    SDL_SetTextureBlendMode(finalTexture, SDL_BLENDMODE_BLEND);
    return finalTexture;
}

void get_component_defaults(ComponentType type, string& value) {
    switch (type) {
        case ComponentType::RESISTOR: value = "1k"; break;
        case ComponentType::CAPACITOR: value = "1u"; break;
        case ComponentType::INDUCTOR: value = "1m"; break;
        case ComponentType::DIODE: value = "1N4148"; break;
        case ComponentType::DC_VOLTAGE_SOURCE: value = "5"; break;
        case ComponentType::DC_CURRENT_SOURCE: value = "1m"; break;
        case ComponentType::AC_VOLTAGE_SOURCE: value = "SIN(0 1 1k)"; break;
        default: value = "?"; break;
    }
}

string generate_component_name(ComponentType type) {
    char prefix = '?';
    switch (type) {
        case ComponentType::RESISTOR: prefix = 'R'; break;
        case ComponentType::CAPACITOR: prefix = 'C'; break;
        case ComponentType::INDUCTOR: prefix = 'L'; break;
        case ComponentType::DIODE: prefix = 'D'; break;
        case ComponentType::DC_VOLTAGE_SOURCE:
        case ComponentType::AC_VOLTAGE_SOURCE: prefix = 'V'; break;
        case ComponentType::DC_CURRENT_SOURCE: prefix = 'I'; break;
        default: break;
    }
    componentNameCounters[type]++;
    return prefix + to_string(componentNameCounters[type]);
}

string component_type_to_string(ComponentType type) {
    switch (type) {
        case ComponentType::RESISTOR: return "RESISTOR";
        case ComponentType::CAPACITOR: return "CAPACITOR";
        case ComponentType::INDUCTOR: return "INDUCTOR";
        case ComponentType::DIODE: return "DIODE";
        case ComponentType::DC_VOLTAGE_SOURCE: return "DC_VOLTAGE_SOURCE";
        case ComponentType::DC_CURRENT_SOURCE: return "DC_CURRENT_SOURCE";
        case ComponentType::AC_VOLTAGE_SOURCE: return "AC_VOLTAGE_SOURCE";
        default: return "UNKNOWN";
    }
}

ComponentType string_to_component_type(const string& str) {
    if (str == "RESISTOR") return ComponentType::RESISTOR;
    if (str == "CAPACITOR") return ComponentType::CAPACITOR;
    if (str == "INDUCTOR") return ComponentType::INDUCTOR;
    if (str == "DIODE") return ComponentType::DIODE;
    if (str == "DC_VOLTAGE_SOURCE") return ComponentType::DC_VOLTAGE_SOURCE;
    if (str == "DC_CURRENT_SOURCE") return ComponentType::DC_CURRENT_SOURCE;
    if (str == "AC_VOLTAGE_SOURCE") return ComponentType::AC_VOLTAGE_SOURCE;
    return ComponentType::NONE;
}

void update_name_counters_from_id(ComponentType type, const string& id) {
    if (id.empty() || !isdigit(id.back())) return;
    size_t first_digit = id.find_first_of("0123456789");
    if (first_digit == string::npos) return;
    int num = stoi(id.substr(first_digit));
    if (num > componentNameCounters[type]) componentNameCounters[type] = num;
}

void save_schematic(const string& filename) {
    ofstream outFile(filename);
    if (!outFile.is_open()) {
        cerr << "Error: Could not open file for saving: " << filename << endl;
        return;
    }
    for (const auto& comp : components) outFile << "COMPONENT," << component_type_to_string(comp.type) << "," << comp.id << "," << comp.node1.x << "," << comp.node1.y << "," << comp.node2.x << "," << comp.node2.y << "," << comp.rotation_angle << "," << comp.value << endl;
    for (const auto& wire : wires) outFile << "WIRE," << wire.start.x << "," << wire.start.y << "," << wire.end.x << "," << wire.end.y << endl;
    for (const auto& label : labels) outFile << "LABEL," << label.text << "," << label.position.x << "," << label.position.y << "," << label.rotation_angle << endl;
    outFile.close();
    currentSchematicFileName = filename;
    schematicModified = false;
    cout << "Schematic saved to " << filename << endl;
}

void load_schematic(const string& filename) {
    ifstream inFile(filename);
    if (!inFile.is_open()) {
        cerr << "Error: Could not open file for loading: " << filename << endl;
        return;
    }
    components.clear();
    wires.clear();
    labels.clear();
    componentNameCounters.clear();
    string line;
    while (getline(inFile, line)) {
        stringstream ss(line);
        string type_str;
        getline(ss, type_str, ',');
        if (type_str == "COMPONENT") {
            string comp_type_str, id_str, x1_str, y1_str, x2_str, y2_str, rot_str, val_str;
            getline(ss, comp_type_str, ',');
            getline(ss, id_str, ',');
            getline(ss, x1_str, ',');
            getline(ss, y1_str, ',');
            getline(ss, x2_str, ',');
            getline(ss, y2_str, ',');
            getline(ss, rot_str, ',');
            getline(ss, val_str);
            ComponentType type = string_to_component_type(comp_type_str);
            components.push_back({ type, { stoi(x1_str), stoi(y1_str) }, { stoi(x2_str), stoi(y2_str) }, id_str, val_str, stoi(rot_str) });
            update_name_counters_from_id(type, id_str);
        }
        else if (type_str == "WIRE") {
            string x1_str, y1_str, x2_str, y2_str;
            getline(ss, x1_str, ',');
            getline(ss, y1_str, ',');
            getline(ss, x2_str, ',');
            getline(ss, y2_str);
            wires.push_back({ { stoi(x1_str), stoi(y1_str) }, { stoi(x2_str), stoi(y2_str) } });
        }
        else if (type_str == "LABEL") {
            string text, x_str, y_str, rot_str;
            getline(ss, text, ',');
            getline(ss, x_str, ',');
            getline(ss, y_str, ',');
            getline(ss, rot_str);
            labels.push_back({ text, { stoi(x_str), stoi(y_str) }, {}, stoi(rot_str) });
        }
    }
    inFile.close();
    currentSchematicFileName = filename;
    schematicModified = false;
}

void draw_schematic_elements(SDL_Renderer* renderer, TTF_Font* valueFont, const vector<ComponentMenuItem>& componentMenuIcons) {
    for (const auto& wire : wires) {
        thickLineRGBA(renderer, wire.start.x, wire.start.y, wire.end.x, wire.end.y, 3, 0xFF, 0xFF, 0xFF, 0xFF);
        filledCircleRGBA(renderer, wire.start.x, wire.start.y, 4, 220, 20, 60, 0xFF);
        filledCircleRGBA(renderer, wire.end.x, wire.end.y, 4, 220, 20, 60, 0xFF);
    }
    for (auto& comp : components) {
        SDL_Texture* iconTexture = nullptr;
        for (const auto& item : componentMenuIcons) if (item.type == comp.type) { iconTexture = item.iconTexture; break; }
        int center_x = (comp.node1.x + comp.node2.x) / 2, center_y = (comp.node1.y + comp.node2.y) / 2;
        if (iconTexture) {
            int icon_size = COMPONENT_DEFAULT_LENGTH;
            SDL_Rect destRect = { center_x - icon_size / 2, center_y - icon_size / 2, icon_size, icon_size };
            SDL_RenderCopyEx(renderer, iconTexture, nullptr, &destRect, comp.rotation_angle, nullptr, SDL_FLIP_NONE);
        }
        filledCircleRGBA(renderer, comp.node1.x, comp.node1.y, 7, 220, 20, 60, 0xFF);
        filledCircleRGBA(renderer, comp.node2.x, comp.node2.y, 7, 220, 20, 60, 0xFF);
        SDL_Color idColor = { 0xAA, 0xAA, 0xFF, 0xFF };
        SDL_Surface* idSurface = TTF_RenderText_Solid(valueFont, comp.id.c_str(), idColor);
        if (idSurface) {
            SDL_Texture* idTexture = SDL_CreateTextureFromSurface(renderer, idSurface);
            comp.id_rect = { center_x - idSurface->w / 2, center_y - 50, idSurface->w, idSurface->h };
            SDL_RenderCopy(renderer, idTexture, nullptr, &comp.id_rect);
            SDL_FreeSurface(idSurface);
            SDL_DestroyTexture(idTexture);
        }
        string text_to_render = (editingComponent == &comp) ? textInputBuffer : comp.value;
        SDL_Color textColor = { 0xFF, 0xFF, 0xFF, 0xFF };
        SDL_Surface* textSurface = TTF_RenderText_Solid(valueFont, text_to_render.c_str(), textColor);
        if (textSurface) {
            SDL_Texture* textTexture = SDL_CreateTextureFromSurface(renderer, textSurface);
            comp.value_rect = { center_x - textSurface->w / 2, center_y + 35, textSurface->w, textSurface->h };
            if (editingComponent == &comp) {
                SDL_SetRenderDrawColor(renderer, 0x00, 0x00, 0x55, 0xFF);
                SDL_RenderFillRect(renderer, &comp.value_rect);
            }
            SDL_RenderCopy(renderer, textTexture, nullptr, &comp.value_rect);
            SDL_FreeSurface(textSurface);
            SDL_DestroyTexture(textTexture);
        }
    }
    for (auto& label : labels) {
        string text_to_render = (editingLabel == &label) ? textInputBuffer : label.text;
        SDL_Color textColor = (label.text == "GND") ? SDL_Color{ 100, 255, 100, 255 } : SDL_Color{ 0x99, 0xFF, 0xFF, 0xFF };
        SDL_Surface* textSurface = TTF_RenderText_Solid(valueFont, text_to_render.c_str(), textColor);
        if (textSurface) {
            SDL_Texture* textTexture = SDL_CreateTextureFromSurface(renderer, textSurface);
            int w = textSurface->w, h = textSurface->h, offsetX = 10, offsetY = 10;
            if (label.rotation_angle == 0) label.text_rect = { label.position.x + offsetX, label.position.y - h - (offsetY / 2), w, h };
            else if (label.rotation_angle == 90) label.text_rect = { label.position.x - w - offsetX, label.position.y - h / 2, w, h };
            else if (label.rotation_angle == 180) label.text_rect = { label.position.x - w / 2, label.position.y + offsetY, w, h };
            else label.text_rect = { label.position.x + offsetX, label.position.y - h / 2, w, h };
            if (editingLabel == &label) {
                SDL_SetRenderDrawColor(renderer, 0x00, 0x33, 0x33, 0xFF);
                SDL_RenderFillRect(renderer, &label.text_rect);
            }
            SDL_RenderCopy(renderer, textTexture, nullptr, &label.text_rect);
            SDL_FreeSurface(textSurface);
            SDL_DestroyTexture(textTexture);
        }
    }
}

// --- New Simulation Handler ---
void run_simulation_in_terminal(DialogType analysis_type, const vector<DialogField>& fields) {
    cout << "\n--- Preparing Simulation ---" << endl;

    // 1. Build SpiceCircuit from GUI data
    auto spice_circuit = make_unique<SpiceEngine::Circuit>();

    // Node Resolution Logic
    map<string, int> point_to_set_id;
    vector<int> parent;
    auto point_to_key = [](SDL_Point p) { return to_string(p.x) + "," + to_string(p.y); };
    set<string> unique_points_keys;
    for (const auto& comp : components) {
        unique_points_keys.insert(point_to_key(comp.node1));
        unique_points_keys.insert(point_to_key(comp.node2));
    }
    for (const auto& wire : wires) {
        unique_points_keys.insert(point_to_key(wire.start));
        unique_points_keys.insert(point_to_key(wire.end));
    }
    for (const auto& label : labels) {
        unique_points_keys.insert(point_to_key(label.position));
    }
    int set_id_counter = 0;
    for (const auto& key : unique_points_keys) {
        point_to_set_id[key] = set_id_counter++;
    }
    parent.resize(set_id_counter);
    iota(parent.begin(), parent.end(), 0);
    function<int(int)> find_set = [&](int i) {
        if (parent[i] == i) return i;
        return parent[i] = find_set(parent[i]);
    };
    auto unite_sets = [&](int a, int b) {
        a = find_set(a); b = find_set(b);
        if (a != b) parent[b] = a;
    };
    for (const auto& wire : wires) {
        unite_sets(point_to_set_id.at(point_to_key(wire.start)), point_to_set_id.at(point_to_key(wire.end)));
    }
    for (const auto& point_key_pair : point_to_set_id) {
        string key = point_key_pair.first;
        stringstream ss(key);
        string x_str, y_str;
        getline(ss, x_str, ',');
        getline(ss, y_str);
        SDL_Point p = { stoi(x_str), stoi(y_str) };
        for (const auto& wire : wires) {
            if (dist_to_segment_sq(p, wire.start, wire.end) < 1.0) {
                unite_sets(point_key_pair.second, point_to_set_id.at(point_to_key(wire.start)));
            }
        }
    }

    // Node Naming Logic
    map<int, string> set_to_node_name;
    int node_name_counter = 1;
    for (const auto& label : labels) {
        if (label.text != "GND") {
            string key = point_to_key(label.position);
            if (point_to_set_id.count(key)) {
                int root = find_set(point_to_set_id.at(key));
                if (!set_to_node_name.count(root)) {
                    set_to_node_name[root] = label.text;
                }
            }
        }
    }
    for (const auto& pair : point_to_set_id) {
        int root = find_set(pair.second);
        if (!set_to_node_name.count(root)) {
            set_to_node_name[root] = "N" + string(3 - to_string(node_name_counter).length(), '0') + to_string(node_name_counter++);
        }
    }
    string ground_node_name = "";
    for (const auto& label : labels) {
        if (label.text == "GND") {
            string key = point_to_key(label.position);
            if (point_to_set_id.count(key)) {
                int root = find_set(point_to_set_id.at(key));
                ground_node_name = set_to_node_name[root];
                spice_circuit->set_ground_node(ground_node_name);
                break;
            }
        }
    }
    auto get_node_name = [&](SDL_Point p) {
        string key = point_to_key(p);
        int root = find_set(point_to_set_id.at(key));
        return set_to_node_name.at(root);
    };

    // Add components to SpiceCircuit
    for (const auto& comp : components) {
        string n1 = get_node_name(comp.node1);
        string n2 = get_node_name(comp.node2);
        switch (comp.type) {
            case ComponentType::RESISTOR: spice_circuit->add_component(make_unique<SpiceEngine::Resistor>(comp.id, n1, n2, comp.value)); break;
            case ComponentType::CAPACITOR: spice_circuit->add_component(make_unique<SpiceEngine::Capacitor>(comp.id, n1, n2, comp.value)); break;
            case ComponentType::INDUCTOR: spice_circuit->add_component(make_unique<SpiceEngine::Inductor>(comp.id, n1, n2, comp.value)); break;
            case ComponentType::DIODE: spice_circuit->add_component(make_unique<SpiceEngine::Diode>(comp.id, n1, n2, comp.value)); break;
            case ComponentType::DC_CURRENT_SOURCE: spice_circuit->add_component(make_unique<SpiceEngine::CurrentSource>(comp.id, n1, n2, comp.value)); break;
            case ComponentType::DC_VOLTAGE_SOURCE:
            case ComponentType::AC_VOLTAGE_SOURCE: {
                vector<string> params;
                string temp_val = comp.value;
                size_t sin_pos = temp_val.find("SIN(");
                if (sin_pos != string::npos) {
                    params.push_back(temp_val); // The parser expects the full string
                } else {
                    params.push_back(comp.value);
                }
                spice_circuit->add_component(make_unique<SpiceEngine::VoltageSource>(comp.id, n1, n2, params));
                break;
            }
            default: break;
        }
    }

    // 2. Run simulation and print results
    try {
        vector<string> vars_to_print;
        if (analysis_type == DialogType::TRANSIENT_ANALYSIS) {
            double t_step = parse_value_with_metric_prefix_util(fields[0].buffer);
            double t_stop = parse_value_with_metric_prefix_util(fields[1].buffer);
            stringstream ss(fields[2].buffer);
            string var;
            while(ss >> var) { vars_to_print.push_back(var); }

            spice_circuit->perform_transient_analysis(t_step, t_stop);

            cout << left << setw(15) << "time (s)";
            for (const auto& v : vars_to_print) cout << setw(15) << v;
            cout << endl;
            for (const auto& result_point : spice_circuit->tran_results) {
                cout << fixed << setprecision(6) << setw(15) << result_point.at("time");
                for (const auto& v : vars_to_print) {
                    try {
                        ostringstream val_ss;
                        val_ss << fixed << setprecision(6) << result_point.at(v);
                        string val_str = val_ss.str();
                        char type = to_lower_util(v)[0];
                        if (type == 'v') val_str += " V";
                        else if (type == 'i') val_str += " A";
                        cout << setw(15) << val_str;
                    } catch (const out_of_range&) {
                        cout << setw(15) << "N/A";
                    }
                }
                cout << endl;
            }

        } else if (analysis_type == DialogType::DC_SWEEP_ANALYSIS) {
            string src_name = fields[0].buffer;
            double start = parse_value_with_metric_prefix_util(fields[1].buffer);
            double end = parse_value_with_metric_prefix_util(fields[2].buffer);
            double inc = parse_value_with_metric_prefix_util(fields[3].buffer);
            stringstream ss(fields[4].buffer);
            string var;
            while(ss >> var) { vars_to_print.push_back(var); }

            spice_circuit->perform_dc_sweep_analysis(src_name, start, end, inc);

            string header_src_name = spice_circuit->dc_sweep_source_name;
            char src_name_char = to_lower_util(header_src_name)[0];
            if (src_name_char == 'v') header_src_name += " (V)";
            else if (src_name_char == 'i') header_src_name += " (A)";
            cout << left << setw(15) << header_src_name;
            for (const auto& v : vars_to_print) cout << setw(15) << v;
            cout << endl;

            for (const auto& result_point : spice_circuit->dc_sweep_results) {
                cout << fixed << setprecision(6) << setw(15) << result_point.at(spice_circuit->dc_sweep_source_name);
                for (const auto& v : vars_to_print) {
                    if (result_point.count(v)) {
                        ostringstream val_ss;
                        val_ss << fixed << setprecision(6) << result_point.at(v);
                        string val_str = val_ss.str();
                        char type = to_lower_util(v)[0];
                        if (type == 'v') val_str += " V";
                        else if (type == 'i') val_str += " A";
                        cout << setw(15) << val_str;
                    } else {
                        cout << setw(15) << "N/A";
                    }
                }
                cout << endl;
            }
        }
    } catch (const exception& e) {
        cerr << "SIMULATION ERROR: " << e.what() << endl;
    }
    cout << "--- Simulation Finished ---" << endl;
}


void close_dialog() {
    activeDialogType = DialogType::NONE;
    dialogFields.clear();
    dialogTitle = "";
    activeDialogFieldIndex = -1;
    currentInteractionMode = InteractionMode::NONE;
    SDL_StopTextInput();
}

void handle_dialog_ok() {
    if (activeDialogType == DialogType::SAVE_AS) {
        if (!dialogFields.empty() && !dialogFields[0].buffer.empty()) save_schematic(SCHEMATICS_PATH + dialogFields[0].buffer);
    }
    else if (activeDialogType == DialogType::OPEN) {
        if (!dialogFields.empty() && !dialogFields[0].buffer.empty()) load_schematic(SCHEMATICS_PATH + dialogFields[0].buffer);
    }
    else if (activeDialogType == DialogType::TRANSIENT_ANALYSIS || activeDialogType == DialogType::DC_SWEEP_ANALYSIS) {
        run_simulation_in_terminal(activeDialogType, dialogFields);
    }
    close_dialog();
}


void setup_dialog(DialogType type, const string& title, const vector<string>& labels) {
    activeDialogType = type;
    dialogTitle = title;
    dialogFields.clear();
    for (const auto& label : labels) dialogFields.push_back({ label, "" });
    activeDialogFieldIndex = 0;
    currentInteractionMode = InteractionMode::DIALOG_ACTIVE;
    SDL_StartTextInput();
}

void render_dialog(SDL_Renderer* renderer, TTF_Font* font) {
    const int dialog_w = 450, dialog_h = 100 + dialogFields.size() * 40 + 50;
    SDL_Rect panelRect = { SCREEN_WIDTH / 2 - dialog_w / 2, SCREEN_HEIGHT / 2 - dialog_h / 2, dialog_w, dialog_h };
    SDL_SetRenderDrawColor(renderer, 0x55, 0x55, 0x55, 0xFF);
    SDL_RenderFillRect(renderer, &panelRect);
    SDL_SetRenderDrawColor(renderer, 0x88, 0x88, 0x88, 0xFF);
    SDL_RenderDrawRect(renderer, &panelRect);
    SDL_Color textColor = { 0xFF, 0xFF, 0xFF, 0xFF };
    SDL_Surface* titleSurface = TTF_RenderText_Solid(font, dialogTitle.c_str(), textColor);
    if (titleSurface) {
        SDL_Texture* titleTexture = SDL_CreateTextureFromSurface(renderer, titleSurface);
        SDL_Rect titleRect = { panelRect.x + (panelRect.w - titleSurface->w) / 2, panelRect.y + 20, titleSurface->w, titleSurface->h };
        SDL_RenderCopy(renderer, titleTexture, nullptr, &titleRect);
        SDL_FreeSurface(titleSurface);
        SDL_DestroyTexture(titleTexture);
    }
    int current_y = panelRect.y + 60;
    for (size_t i = 0; i < dialogFields.size(); ++i) {
        SDL_Surface* labelSurface = TTF_RenderText_Solid(font, dialogFields[i].label.c_str(), textColor);
        if (labelSurface) {
            SDL_Texture* labelTexture = SDL_CreateTextureFromSurface(renderer, labelSurface);
            SDL_Rect labelRect = { panelRect.x + 20, current_y + (30 - labelSurface->h) / 2, labelSurface->w, labelSurface->h };
            SDL_RenderCopy(renderer, labelTexture, nullptr, &labelRect);
            SDL_FreeSurface(labelSurface);
            SDL_DestroyTexture(labelTexture);
        }
        SDL_Rect inputRect = { panelRect.x + 150, current_y, panelRect.w - 170, 30 };
        dialogFields[i].input_rect = inputRect;
        SDL_SetRenderDrawColor(renderer, (activeDialogFieldIndex == i) ? 0x00 : 0x33, (activeDialogFieldIndex == i) ? 0x00 : 0x33, (activeDialogFieldIndex == i) ? 0x44 : 0x33, 0xFF);
        SDL_RenderFillRect(renderer, &inputRect);
        SDL_SetRenderDrawColor(renderer, 0xAAAAAA, 0xAAAAAA, 0xAAAAAA, 0xFF);
        SDL_RenderDrawRect(renderer, &inputRect);
        if (!dialogFields[i].buffer.empty()) {
            SDL_Surface* textSurface = TTF_RenderText_Solid(font, dialogFields[i].buffer.c_str(), textColor);
            if (textSurface) {
                SDL_Texture* textTexture = SDL_CreateTextureFromSurface(renderer, textSurface);
                SDL_Rect textRect = { inputRect.x + 5, inputRect.y + (inputRect.h - textSurface->h) / 2, textSurface->w, textSurface->h };
                SDL_RenderCopy(renderer, textTexture, nullptr, &textRect);
                SDL_FreeSurface(textSurface);
                SDL_DestroyTexture(textTexture);
            }
        }
        current_y += 40;
    }
    Button okButton("OK", panelRect.x + panelRect.w - 190, panelRect.y + panelRect.h - 45, 80, 30, handle_dialog_ok);
    Button cancelButton("Cancel", panelRect.x + panelRect.w - 100, panelRect.y + panelRect.h - 45, 80, 30, close_dialog);
    int mx, my;
    SDL_GetMouseState(&mx, &my);
    okButton.m_is_hovered = SDL_PointInRect(new const SDL_Point{ mx, my }, &okButton.m_position);
    cancelButton.m_is_hovered = SDL_PointInRect(new const SDL_Point{ mx, my }, &cancelButton.m_position);
    okButton.render(renderer, font);
    cancelButton.render(renderer, font);
}

// --- Main Function ---
int main(int argc, char* args[]) {
    if (SDL_Init(SDL_INIT_VIDEO) < 0 || TTF_Init() == -1) {
        cerr << "Core SDL init failed." << endl;
        return -1;
    }
    if (!(IMG_Init(IMG_INIT_PNG) & IMG_INIT_PNG)) {
        cerr << "SDL_image init failed." << endl;
        return -1;
    }
    SDL_Window* window = SDL_CreateWindow("SUTSpice | Integrated Simulator", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, SCREEN_WIDTH, SCREEN_HEIGHT, SDL_WINDOW_SHOWN);
    SDL_Renderer* renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED | SDL_RENDERER_PRESENTVSYNC);
    TTF_Font* uiFont = TTF_OpenFont("C:\\Windows\\Fonts\\arial.ttf", 16);
    TTF_Font* valueFont = TTF_OpenFont("C:\\Windows\\Fonts\\arialbd.ttf", 18);
    if (!window || !renderer || !uiFont || !valueFont) {
        cerr << "Window or font creation failed." << endl;
        return -1;
    }

    // ... (The rest of the main function remains largely the same, setting up buttons and the main loop)
    vector<Button> topBarButtons;
    topBarButtons.emplace_back("File", 10, 5, 80, 30, []() { isFileMenuOpen = !isFileMenuOpen; isComponentsMenuOpen = false; isSimulateMenuOpen = false; currentInteractionMode = InteractionMode::NONE; });
    topBarButtons.emplace_back("Simulate", 100, 5, 100, 30, []() { isSimulateMenuOpen = !isSimulateMenuOpen; isFileMenuOpen = false; isComponentsMenuOpen = false; currentInteractionMode = InteractionMode::NONE; });
    topBarButtons.emplace_back("Add Component", 210, 5, 150, 30, []() { isComponentsMenuOpen = !isComponentsMenuOpen; isFileMenuOpen = false; isSimulateMenuOpen = false; currentInteractionMode = InteractionMode::NONE; });
    topBarButtons.emplace_back("Add Wire", 370, 5, 100, 30, []() { currentInteractionMode = InteractionMode::WIRING; });
    topBarButtons.emplace_back("Add GND", 480, 5, 100, 30, []() { currentInteractionMode = InteractionMode::PLACE_GND_LABEL; });
    topBarButtons.emplace_back("Add Label", 590, 5, 100, 30, []() { currentInteractionMode = InteractionMode::PLACE_LABEL; });
    topBarButtons.emplace_back("Delete", 700, 5, 100, 30, []() { currentInteractionMode = InteractionMode::DELETE_ITEM; });

    vector<Button> fileMenuButtons;
    fileMenuButtons.emplace_back("New Schematic", 10, 40, 180, 30, []() { components.clear(); wires.clear(); labels.clear(); componentNameCounters.clear(); currentSchematicFileName = ""; schematicModified = false; isFileMenuOpen = false; });
    fileMenuButtons.emplace_back("Open...", 10, 70, 180, 30, []() { setup_dialog(DialogType::OPEN, "Open Schematic", { "Filename" }); isFileMenuOpen = false; });
    fileMenuButtons.emplace_back("Save", 10, 100, 180, 30, []() { if (currentSchematicFileName.empty()) { setup_dialog(DialogType::SAVE_AS, "Save Schematic As", { "Filename" }); } else { save_schematic(currentSchematicFileName); } isFileMenuOpen = false; });
    fileMenuButtons.emplace_back("Save As...", 10, 130, 180, 30, []() { setup_dialog(DialogType::SAVE_AS, "Save Schematic As", { "Filename" }); isFileMenuOpen = false; });
    bool quit_flag = false;
    fileMenuButtons.emplace_back("Exit", 10, 160, 180, 30, [&quit_flag]() { quit_flag = true; });

    vector<Button> simulateMenuButtons;
    simulateMenuButtons.emplace_back("DC Sweep Analysis", 100, 40, 180, 30, []() { setup_dialog(DialogType::DC_SWEEP_ANALYSIS, "DC Sweep Analysis", { "Source Name", "Start Value", "End Value", "Increment", "Wanted Values" }); isSimulateMenuOpen = false; });
    simulateMenuButtons.emplace_back("Transient Analysis", 100, 70, 180, 30, []() { setup_dialog(DialogType::TRANSIENT_ANALYSIS, "Transient Analysis", { "Tstep", "Tstop", "Wanted Values" }); isSimulateMenuOpen = false; });

    vector<ComponentMenuItem> componentMenuItems;
    const int COMP_ITEM_W = 150, COMP_ITEM_H = 100, COMP_ITEM_PAD = 10;
    auto selectComp = [](ComponentType t) { selectedComponentType = t; currentInteractionMode = InteractionMode::PLACE_COMPONENT; isComponentsMenuOpen = false; placementRotation = 0; };
    componentMenuItems.push_back({ "Resistor", ComponentType::RESISTOR, loadAndProcessTexture(ASSET_PATH + "resistor.png", renderer), Button("", 210, 40, COMP_ITEM_W, COMP_ITEM_H, [=]() { selectComp(ComponentType::RESISTOR); }) });
    componentMenuItems.push_back({ "Capacitor", ComponentType::CAPACITOR, loadAndProcessTexture(ASSET_PATH + "capacitor.png", renderer), Button("", 210 + (COMP_ITEM_W + COMP_ITEM_PAD), 40, COMP_ITEM_W, COMP_ITEM_H, [=]() { selectComp(ComponentType::CAPACITOR); }) });
    componentMenuItems.push_back({ "Inductor", ComponentType::INDUCTOR, loadAndProcessTexture(ASSET_PATH + "inductor.png", renderer), Button("", 210 + 2 * (COMP_ITEM_W + COMP_ITEM_PAD), 40, COMP_ITEM_W, COMP_ITEM_H, [=]() { selectComp(ComponentType::INDUCTOR); }) });
    componentMenuItems.push_back({ "Diode", ComponentType::DIODE, loadAndProcessTexture(ASSET_PATH + "diode.png", renderer), Button("", 210 + 3 * (COMP_ITEM_W + COMP_ITEM_PAD), 40, COMP_ITEM_W, COMP_ITEM_H, [=]() { selectComp(ComponentType::DIODE); }) });
    componentMenuItems.push_back({ "DC Voltage Src", ComponentType::DC_VOLTAGE_SOURCE, loadAndProcessTexture(ASSET_PATH + "voltage_source.png", renderer), Button("", 210, 40 + (COMP_ITEM_H + COMP_ITEM_PAD), COMP_ITEM_W, COMP_ITEM_H, [=]() { selectComp(ComponentType::DC_VOLTAGE_SOURCE); }) });
    componentMenuItems.push_back({ "DC Current Src", ComponentType::DC_CURRENT_SOURCE, loadAndProcessTexture(ASSET_PATH + "current_source.png", renderer), Button("", 210 + (COMP_ITEM_W + COMP_ITEM_PAD), 40 + (COMP_ITEM_H + COMP_ITEM_PAD), COMP_ITEM_W, COMP_ITEM_H, [=]() { selectComp(ComponentType::DC_CURRENT_SOURCE); }) });
    componentMenuItems.push_back({ "AC Voltage Src", ComponentType::AC_VOLTAGE_SOURCE, loadAndProcessTexture(ASSET_PATH + "ac_voltage_source.png", renderer), Button("", 210 + 2 * (COMP_ITEM_W + COMP_ITEM_PAD), 40 + (COMP_ITEM_H + COMP_ITEM_PAD), COMP_ITEM_W, COMP_ITEM_H, [=]() { selectComp(ComponentType::AC_VOLTAGE_SOURCE); }) });

    SDL_Event e;
    while (!quit_flag) {
        // Main event loop
        while (SDL_PollEvent(&e) != 0) {
            if (e.type == SDL_QUIT) quit_flag = true;
            // ... (rest of the event handling logic)
            if (currentInteractionMode == InteractionMode::DIALOG_ACTIVE) {
                if (e.type == SDL_KEYDOWN) {
                    if (e.key.keysym.sym == SDLK_RETURN) handle_dialog_ok();
                    else if (e.key.keysym.sym == SDLK_ESCAPE) close_dialog();
                    else if (e.key.keysym.sym == SDLK_TAB && activeDialogFieldIndex != -1) activeDialogFieldIndex = (activeDialogFieldIndex + 1) % dialogFields.size();
                    else if (e.key.keysym.sym == SDLK_BACKSPACE && activeDialogFieldIndex != -1 && !dialogFields[activeDialogFieldIndex].buffer.empty()) dialogFields[activeDialogFieldIndex].buffer.pop_back();
                } else if (e.type == SDL_TEXTINPUT && activeDialogFieldIndex != -1) dialogFields[activeDialogFieldIndex].buffer += e.text.text;
                else if (e.type == SDL_MOUSEBUTTONDOWN) {
                    int mx, my; SDL_GetMouseState(&mx, &my);
                    for (size_t i = 0; i < dialogFields.size(); ++i) if (SDL_PointInRect(new const SDL_Point{ mx, my }, &dialogFields[i].input_rect)) { activeDialogFieldIndex = i; break; }
                    const int dialog_w = 450, dialog_h = 100 + dialogFields.size() * 40 + 50;
                    SDL_Rect panelRect = { SCREEN_WIDTH / 2 - dialog_w / 2, SCREEN_HEIGHT / 2 - dialog_h / 2, dialog_w, dialog_h };
                    SDL_Rect okRect = { panelRect.x + panelRect.w - 190, panelRect.y + panelRect.h - 45, 80, 30 };
                    SDL_Rect cancelRect = { panelRect.x + panelRect.w - 100, panelRect.y + panelRect.h - 45, 80, 30 };
                    if (SDL_PointInRect(new const SDL_Point{ mx, my }, &okRect)) handle_dialog_ok();
                    else if (SDL_PointInRect(new const SDL_Point{ mx, my }, &cancelRect)) close_dialog();
                }
                continue;
            }
            if (currentInteractionMode == InteractionMode::EDITING_COMPONENT_VALUE || currentInteractionMode == InteractionMode::EDITING_LABEL_TEXT) {
                if (e.type == SDL_KEYDOWN) {
                    if (e.key.keysym.sym == SDLK_RETURN || e.key.keysym.sym == SDLK_ESCAPE) {
                        if (e.key.keysym.sym == SDLK_RETURN) {
                            if (editingComponent) editingComponent->value = textInputBuffer;
                            if (editingLabel) editingLabel->text = textInputBuffer;
                        }
                        else if (editingLabel && editingLabel->text.empty()) labels.pop_back();
                        editingComponent = nullptr; editingLabel = nullptr;
                        currentInteractionMode = InteractionMode::NONE;
                        SDL_StopTextInput();
                    }
                    else if (e.key.keysym.sym == SDLK_BACKSPACE && !textInputBuffer.empty()) textInputBuffer.pop_back();
                    else if (e.key.keysym.sym == SDLK_r && (e.key.keysym.mod & KMOD_CTRL) && editingLabel) editingLabel->rotation_angle = (editingLabel->rotation_angle + 90) % 360;
                }
                else if (e.type == SDL_TEXTINPUT) textInputBuffer += e.text.text;
                continue;
            }
            if (e.type == SDL_KEYDOWN) {
                if (e.key.keysym.sym == SDLK_ESCAPE) { currentInteractionMode = InteractionMode::NONE; isDrawingWire = false; }
                else if (e.key.keysym.sym == SDLK_r && currentInteractionMode == InteractionMode::PLACE_COMPONENT) placementRotation = (placementRotation + 90) % 360;
            }
            if (e.type == SDL_MOUSEBUTTONDOWN) {
                int mouseX, mouseY; SDL_GetMouseState(&mouseX, &mouseY);
                SDL_Point snappedPoint = snap_to_grid(mouseX, mouseY);
                if (e.button.button == SDL_BUTTON_LEFT) {
                    bool itemClicked = false;
                    for (auto& comp : components) if (SDL_PointInRect(&snappedPoint, &comp.value_rect)) { editingComponent = &comp; textInputBuffer = comp.value; currentInteractionMode = InteractionMode::EDITING_COMPONENT_VALUE; SDL_StartTextInput(); itemClicked = true; break; }
                    if (itemClicked) continue;
                    if (currentInteractionMode == InteractionMode::PLACE_LABEL) { labels.push_back({ "", snappedPoint }); editingLabel = &labels.back(); textInputBuffer = ""; currentInteractionMode = InteractionMode::EDITING_LABEL_TEXT; SDL_StartTextInput(); }
                    else if (currentInteractionMode == InteractionMode::PLACE_GND_LABEL) { labels.push_back({ "GND", snappedPoint, {}, 0 }); schematicModified = true; currentInteractionMode = InteractionMode::NONE; }
                    else if (currentInteractionMode == InteractionMode::WIRING) { firstInteractionPoint = snappedPoint; isDrawingWire = true; }
                    else if (currentInteractionMode == InteractionMode::PLACE_COMPONENT) {
                        Component newComp = { selectedComponentType, {}, {}, generate_component_name(selectedComponentType), "", placementRotation };
                        get_component_defaults(newComp.type, newComp.value);
                        int half_len = COMPONENT_DEFAULT_LENGTH / 2;
                        if (placementRotation == 0) { newComp.node1 = { snappedPoint.x - half_len, snappedPoint.y }; newComp.node2 = { snappedPoint.x + half_len, snappedPoint.y }; }
                        else if (placementRotation == 90) { newComp.node1 = { snappedPoint.x, snappedPoint.y - half_len }; newComp.node2 = { snappedPoint.x, snappedPoint.y + half_len }; }
                        else if (placementRotation == 180) { newComp.node1 = { snappedPoint.x + half_len, snappedPoint.y }; newComp.node2 = { snappedPoint.x - half_len, snappedPoint.y }; }
                        else { newComp.node1 = { snappedPoint.x, snappedPoint.y + half_len }; newComp.node2 = { snappedPoint.x, snappedPoint.y - half_len }; }
                        components.push_back(newComp); schematicModified = true;
                    }
                    else if (currentInteractionMode == InteractionMode::DELETE_ITEM) {
                        const double DELETE_THRESHOLD_SQ = 10 * 10;
                        SDL_Point clickPoint = { mouseX, mouseY };
                        int wire_to_delete = -1, component_to_delete = -1, label_to_delete = -1;
                        double min_dist_sq = DELETE_THRESHOLD_SQ;
                        for (int i = 0; i < labels.size(); ++i) if (SDL_PointInRect(&clickPoint, &labels[i].text_rect)) { label_to_delete = i; break; }
                        if (label_to_delete == -1) {
                            for (int i = 0; i < wires.size(); ++i) { double d = dist_to_segment_sq(clickPoint, wires[i].start, wires[i].end); if (d < min_dist_sq) { min_dist_sq = d; wire_to_delete = i; component_to_delete = -1; } }
                            for (int i = 0; i < components.size(); ++i) { double d = dist_to_segment_sq(clickPoint, components[i].node1, components[i].node2); if (d < min_dist_sq) { min_dist_sq = d; component_to_delete = i; wire_to_delete = -1; } }
                        }
                        if (label_to_delete != -1) { labels.erase(labels.begin() + label_to_delete); schematicModified = true; }
                        else if (wire_to_delete != -1) { wires.erase(wires.begin() + wire_to_delete); schematicModified = true; }
                        else if (component_to_delete != -1) { components.erase(components.begin() + component_to_delete); schematicModified = true; }
                    }
                }
                else if (e.button.button == SDL_BUTTON_RIGHT) { currentInteractionMode = InteractionMode::NONE; isDrawingWire = false; }
            }
            else if (e.type == SDL_MOUSEBUTTONUP && e.button.button == SDL_BUTTON_LEFT && currentInteractionMode == InteractionMode::WIRING && isDrawingWire) {
                int mouseX, mouseY; SDL_GetMouseState(&mouseX, &mouseY);
                SDL_Point snappedPoint = snap_to_grid(mouseX, mouseY);
                if (snappedPoint.x != firstInteractionPoint.x || snappedPoint.y != firstInteractionPoint.y) { wires.push_back({ firstInteractionPoint, snappedPoint }); schematicModified = true; }
                isDrawingWire = false;
            }
            if (isFileMenuOpen) for (auto& b : fileMenuButtons) b.handle_event(&e);
            if (isComponentsMenuOpen) for (auto& i : componentMenuItems) i.button.handle_event(&e);
            if (isSimulateMenuOpen) for (auto& b : simulateMenuButtons) b.handle_event(&e);
            for (auto& b : topBarButtons) b.handle_event(&e);
        }

        // --- Rendering ---
        SDL_SetRenderDrawColor(renderer, 0x22, 0x22, 0x22, 0xFF);
        SDL_RenderClear(renderer);
        draw_grid(renderer);
        draw_schematic_elements(renderer, valueFont, componentMenuItems);

        if (currentInteractionMode == InteractionMode::WIRING && isDrawingWire) {
            int mx, my; SDL_GetMouseState(&mx, &my);
            thickLineRGBA(renderer, firstInteractionPoint.x, firstInteractionPoint.y, snap_to_grid(mx, my).x, snap_to_grid(mx, my).y, 3, 0xFF, 0xFF, 0xFF, 0xFF);
        }
        if (currentInteractionMode == InteractionMode::PLACE_COMPONENT) {
            int mx, my; SDL_GetMouseState(&mx, &my);
            SDL_Point center = snap_to_grid(mx, my);
            SDL_Texture* icon = nullptr;
            for (const auto& item : componentMenuItems) if (item.type == selectedComponentType) icon = item.iconTexture;
            if (icon) {
                SDL_SetTextureAlphaMod(icon, 150);
                SDL_Rect dest = { center.x - COMPONENT_DEFAULT_LENGTH / 2, center.y - COMPONENT_DEFAULT_LENGTH / 2, COMPONENT_DEFAULT_LENGTH, COMPONENT_DEFAULT_LENGTH };
                SDL_RenderCopyEx(renderer, icon, nullptr, &dest, placementRotation, nullptr, SDL_FLIP_NONE);
                SDL_SetTextureAlphaMod(icon, 255);
            }
        }

        SDL_Rect topBarRect = { 0, 0, SCREEN_WIDTH, TOP_BAR_HEIGHT };
        SDL_SetRenderDrawColor(renderer, 0x33, 0x33, 0x33, 0xFF);
        SDL_RenderFillRect(renderer, &topBarRect);
        for (auto& b : topBarButtons) b.render(renderer, uiFont);

        if (isFileMenuOpen) for (auto& b : fileMenuButtons) b.render(renderer, uiFont);
        if (isSimulateMenuOpen) for (auto& b : simulateMenuButtons) b.render(renderer, uiFont);
        if (isComponentsMenuOpen) {
            SDL_Rect menuPanel = { 205, 35, 4 * (COMP_ITEM_W + COMP_ITEM_PAD) + 5, 2 * (COMP_ITEM_H + COMP_ITEM_PAD) + 5 };
            SDL_SetRenderDrawColor(renderer, 0x3A, 0x3A, 0x3A, 0xFF);
            SDL_RenderFillRect(renderer, &menuPanel);
            for (auto& item : componentMenuItems) {
                item.button.render(renderer, uiFont);
                if (item.iconTexture) {
                    SDL_Rect iconRect = { item.button.m_position.x + 45, item.button.m_position.y + 10, 60, 60 };
                    SDL_RenderCopy(renderer, item.iconTexture, nullptr, &iconRect);
                }
                SDL_Surface* textSurface = TTF_RenderText_Solid(uiFont, item.name.c_str(), { 0xFF, 0xFF, 0xFF, 0xFF });
                if (textSurface) {
                    SDL_Texture* textTexture = SDL_CreateTextureFromSurface(renderer, textSurface);
                    SDL_Rect textRect = { item.button.m_position.x + (item.button.m_position.w - textSurface->w) / 2, item.button.m_position.y + 75, textSurface->w, textSurface->h };
                    SDL_RenderCopy(renderer, textTexture, nullptr, &textRect);
                    SDL_FreeSurface(textSurface);
                    SDL_DestroyTexture(textTexture);
                }
            }
        }
        if (currentInteractionMode == InteractionMode::DIALOG_ACTIVE) render_dialog(renderer, uiFont);
        SDL_RenderPresent(renderer);
    }

    // --- Cleanup ---
    SDL_StopTextInput();
    for (auto& item : componentMenuItems) if (item.iconTexture) SDL_DestroyTexture(item.iconTexture);
    TTF_CloseFont(uiFont);
    TTF_CloseFont(valueFont);
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    IMG_Quit();
    TTF_Quit();
    SDL_Quit();
    return 0;
}
