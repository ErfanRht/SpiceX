#include <SDL2/SDL.h>
#include <SDL2/SDL_image.h>
#include <SDL2/SDL_ttf.h>
#ifdef _WIN32
#include <SDL2/SDL2_gfx.h>
#include <windows.h>
#else
#include <SDL2/SDL2_gfxPrimitives.h>
#include <dirent.h>
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
#include <cfloat>
#include <complex>
using namespace std;
const int SCREEN_WIDTH = 1280;
const int SCREEN_HEIGHT = 720;
const int TOP_BAR_HEIGHT = 40;
const int GRID_SPACING = 20;
const int COMPONENT_DEFAULT_LENGTH = GRID_SPACING * 4;
const double PI = 3.14159265358979323846;
const string BASE_PATH = "C:/Users/Erfan/Dev/Cpp/sutSpice_phase2/";
const string ASSET_PATH = BASE_PATH + "assets/";
const string SCHEMATICS_PATH = BASE_PATH + "schematics/";
enum class AppState {
    SCHEMATIC_EDITOR,
    RESULTS_VIEW
};
AppState currentAppState = AppState::SCHEMATIC_EDITOR;
struct PlottedVariable {
    string name;
    SDL_Color color;
};
vector<PlottedVariable> plotted_variables;
struct ScaleSettings {
    string x_min_str, x_max_str, y_min_str, y_max_str;
};
ScaleSettings manual_scale;
bool show_scale_dialog = false;
int cursor_index = 0;
map<string, pair<string, string>> component_id_to_nodes;
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
string format_value_with_metric_prefix(double value) {
    if (abs(value) < 1e-15) {
        return "0.0";
    }
    const char* prefixes[] = {"p", "n", "u", "m", "", "k", "M", "G"};
    int prefix_index = 4;
    double abs_val = abs(value);
    if (abs_val >= 1.0) {
        while (abs_val >= 1000.0 && prefix_index < 7) {
            abs_val /= 1000.0;
            prefix_index++;
        }
    } else {
        while (abs_val < 1.0 && prefix_index > 0) {
            abs_val *= 1000.0;
            prefix_index--;
        }
    }
    ostringstream oss;
    oss << fixed << setprecision(2) << (value < 0 ? -abs_val : abs_val) << prefixes[prefix_index];
    return oss.str();
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
vector<complex<double>> complex_gaussian_elimination(vector<vector<complex<double>>> A, vector<complex<double>> b) {
    int n = A.size();
    if (n == 0 || (n > 0 && (A[0].size() != n || b.size() != n))) {
        throw runtime_error("Invalid complex matrix or vector dimensions.");
    }
    for (int i = 0; i < n; ++i) {
        int max_row = i;
        for (int k = i + 1; k < n; ++k) {
            if (abs(A[k][i]) > abs(A[max_row][i])) {
                max_row = k;
            }
        }
        swap(A[i], A[max_row]);
        swap(b[i], b[max_row]);
        if (abs(A[i][i]) < 1e-12) continue;
        for (int k = i + 1; k < n; ++k) {
            complex<double> factor = A[k][i] / A[i][i];
            for (int j = i; j < n; ++j) {
                A[k][j] -= factor * A[i][j];
            }
            b[k] -= factor * b[i];
        }
    }
    vector<complex<double>> x(n);
    for (int i = n - 1; i >= 0; --i) {
        if (abs(A[i][i]) < 1e-12) {
            if (abs(b[i]) > 1e-12) {
                throw runtime_error("Complex system is inconsistent.");
            }
            x[i] = 0;
        } else {
            x[i] = b[i];
            for (int j = i + 1; j < n; ++j) {
                x[i] -= A[i][j] * x[j];
            }
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
vector<string> get_schematic_files() {
    vector<string> files;
    string path = SCHEMATICS_PATH;
#ifdef _WIN32
    string search_path = path + "*.txt";
    WIN32_FIND_DATAA find_data;
    HANDLE h_find = FindFirstFileA(search_path.c_str(), &find_data);
    if (h_find == INVALID_HANDLE_VALUE) {
        files.push_back("Error: Directory not found or empty.");
        return files;
    }
    do {
        if (!(find_data.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY)) {
            files.push_back(find_data.cFileName);
        }
    } while (FindNextFileA(h_find, &find_data) != 0);
    FindClose(h_find);
#else
    DIR *dir;
    struct dirent *ent;
    if ((dir = opendir(path.c_str())) != NULL) {
        while ((ent = readdir(dir)) != NULL) {
            string filename = ent->d_name;
            if (filename.length() > 4 && filename.substr(filename.length() - 4) == ".txt") {
                 files.push_back(filename);
            }
        }
        closedir(dir);
    } else {
        files.push_back("Error: Could not read directory.");
    }
#endif
    if (files.empty()) {
        files.push_back("No .txt files found.");
    }
    return files;
}
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
        virtual void ac_stamp(Circuit&, vector<vector<complex<double>>>&, vector<vector<complex<double>>>&, vector<vector<complex<double>>>&, vector<vector<complex<double>>>&, vector<complex<double>>&, vector<complex<double>>&, map<string, int>&, double) = 0;
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
        void ac_stamp(Circuit& circuit, vector<vector<complex<double>>>& G, vector<vector<complex<double>>>&, vector<vector<complex<double>>>&, vector<vector<complex<double>>>&, vector<complex<double>>&, vector<complex<double>>&, map<string, int>&, double) override;
    };
    class VoltageSource : public Component {
    public:
        enum SourceType { DC, SINUSOIDAL, PULSE, AC };
        SourceType sourceType;
        double dc_offset, amplitude, frequency, phase_degrees;
        double v1, v2, td, tr, tf, pw, per;
        vector<string> raw_params;
        VoltageSource(const string& v_name, const string& n1, const string& n2, const vector<string>& params);
        ComponentType get_type() const override { return ComponentType::VoltageSource; }
        string to_netlist_string() const override;
        void stamp(Circuit& circuit, vector<vector<double>>&, vector<vector<double>>& B, vector<vector<double>>& C, vector<vector<double>>&, vector<double>&, vector<double>& E, map<string, int>& m_map, double, const vector<double>&) override;
        void ac_stamp(Circuit& circuit, vector<vector<complex<double>>>&, vector<vector<complex<double>>>& B, vector<vector<complex<double>>>& C, vector<vector<complex<double>>>&, vector<complex<double>>&, vector<complex<double>>& E, map<string, int>& m_map, double) override;
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
        void ac_stamp(Circuit& circuit, vector<vector<complex<double>>>& G, vector<vector<complex<double>>>&, vector<vector<complex<double>>>&, vector<vector<complex<double>>>&, vector<complex<double>>&, vector<complex<double>>&, map<string, int>&, double omega) override;
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
        void ac_stamp(Circuit& circuit, vector<vector<complex<double>>>&, vector<vector<complex<double>>>& B, vector<vector<complex<double>>>& C, vector<vector<complex<double>>>& D, vector<complex<double>>&, vector<complex<double>>& E, map<string, int>& m_map, double omega) override;
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
        void ac_stamp(Circuit& circuit, vector<vector<complex<double>>>&, vector<vector<complex<double>>>&, vector<vector<complex<double>>>&, vector<vector<complex<double>>>&, vector<complex<double>>& J, vector<complex<double>>&, map<string, int>&, double) override;
    };
    class Diode : public Component {
    public:
        bool is_on = false;
        const double Ron = 1e-3;
        const double Roff = 1e9;
        Diode(const string& name, const string& n1, const string& n2, const string& model) : Component(name, n1, n2, 0.0) {
            if (to_lower_util(model) != "ideal" && to_lower_util(model) != "1n4148") {
                throw runtime_error("Only 'ideal' or '1N4148' diode model is supported. Got: " + model);
            }
        }
        ComponentType get_type() const override { return ComponentType::Diode; }
        string to_netlist_string() const override { return name + " " + node1_name + " " + node2_name + " ideal"; }
        void stamp(Circuit& circuit, vector<vector<double>>& G, vector<vector<double>>&, vector<vector<double>>&, vector<vector<double>>&, vector<double>&, vector<double>&, map<string, int>&, double, const vector<double>&) override;
        void ac_stamp(Circuit& circuit, vector<vector<complex<double>>>& G, vector<vector<complex<double>>>&, vector<vector<complex<double>>>&, vector<vector<complex<double>>>&, vector<complex<double>>&, vector<complex<double>>&, map<string, int>&, double) override;
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
        double tran_t_stop = 0.0;
        bool dc_sweep_solved = false;
        vector<ResultPoint> dc_sweep_results;
        string dc_sweep_source_name;
        double dc_start = 0.0, dc_end = 0.0;
        bool ac_sweep_solved = false;
        vector<ResultPoint> ac_sweep_results;
        bool phase_sweep_solved = false;
        vector<ResultPoint> phase_sweep_results;
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
        void build_ac_mna_matrix(vector<vector<complex<double>>>& A, vector<complex<double>>& z, double omega);
        void calculate_and_store_passive_currents(ResultPoint& result_point, const ResultPoint& prev_result_point, double h);
        void perform_transient_analysis(double t_step, double t_stop);
        void perform_dc_sweep_analysis(const string& src, double start, double end, double inc);
        void perform_ac_sweep_analysis(const string& src_name, double f_start, double f_end, int points_per_decade);
        void perform_phase_sweep_analysis(double freq, const string& src, double p_start, double p_end, double p_inc);
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
    void Resistor::ac_stamp(Circuit& circuit, vector<vector<complex<double>>>& G, vector<vector<complex<double>>>&, vector<vector<complex<double>>>&, vector<vector<complex<double>>>&, vector<complex<double>>&, vector<complex<double>>&, map<string, int>&, double) {
        complex<double> conductance = 1.0 / value;
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
            : Component(v_name, n1, n2, 0.0), sourceType(DC), dc_offset(0.0), amplitude(0.0), frequency(0.0), phase_degrees(0.0), v1(0), v2(0), td(0), tr(0), tf(0), pw(0), per(0), raw_params(params) {
        if (params.empty()) throw runtime_error("No value or parameters provided for voltage source " + name);
        string first_param_upper = to_lower_util(params[0]);
        if (first_param_upper.rfind("sin", 0) == 0) {
            string combined_params = params[0];
            for (size_t i = 1; i < params.size(); ++i) {
                combined_params += " " + params[i];
            }
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
            phase_degrees = 0.0;
            this->value = dc_offset;
        }
        else if (first_param_upper.rfind("pulse", 0) == 0) {
            string combined_params = params[0];
            for (size_t i = 1; i < params.size(); ++i) {
                combined_params += " " + params[i];
            }
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
        else if (first_param_upper == "ac") {
            sourceType = AC;
            if (params.size() < 2) throw runtime_error("AC source requires amplitude and optional phase.");
            amplitude = parse_value_with_metric_prefix_util(params[1]);
            phase_degrees = (params.size() > 2) ? parse_value_with_metric_prefix_util(params[2]) : 0.0;
            this->value = 0;
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
    void VoltageSource::ac_stamp(Circuit& circuit, vector<vector<complex<double>>>&, vector<vector<complex<double>>>& B, vector<vector<complex<double>>>& C, vector<vector<complex<double>>>&, vector<complex<double>>&, vector<complex<double>>& E, map<string, int>& m_map, double) {
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
        if (sourceType == AC || sourceType == SINUSOIDAL) {
            E[m_idx] += polar(amplitude, phase_degrees * PI / 180.0);
        }
    }
    void VoltageSource::update_time_dependant_value(double time) {
        if (sourceType == SINUSOIDAL) {
            value = (frequency > 0) ? (dc_offset + amplitude * sin(2 * PI * frequency * time)) : (dc_offset + amplitude);
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
    void CurrentSource::ac_stamp(Circuit& circuit, vector<vector<complex<double>>>&, vector<vector<complex<double>>>&, vector<vector<complex<double>>>&, vector<vector<complex<double>>>&, vector<complex<double>>& J, vector<complex<double>>&, map<string, int>&, double) {
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
    void Capacitor::ac_stamp(Circuit& circuit, vector<vector<complex<double>>>& G, vector<vector<complex<double>>>&, vector<vector<complex<double>>>&, vector<vector<complex<double>>>&, vector<complex<double>>&, vector<complex<double>>&, map<string, int>&, double omega) {
        complex<double> admittance(0, omega * value);
        int idx1 = circuit.get_node_matrix_index(node1_name);
        int idx2 = circuit.get_node_matrix_index(node2_name);
        if (idx1 >= 0) G[idx1][idx1] += admittance;
        if (idx2 >= 0) G[idx2][idx2] += admittance;
        if (idx1 >= 0 && idx2 >= 0) {
            G[idx1][idx2] -= admittance;
            G[idx2][idx1] -= admittance;
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
    void Inductor::ac_stamp(Circuit& circuit, vector<vector<complex<double>>>&, vector<vector<complex<double>>>& B, vector<vector<complex<double>>>& C, vector<vector<complex<double>>>& D, vector<complex<double>>&, vector<complex<double>>&, map<string, int>& m_map, double omega) {
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
        D[m_idx][m_idx] -= complex<double>(0, omega * value);
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
    void Diode::ac_stamp(Circuit& circuit, vector<vector<complex<double>>>& G, vector<vector<complex<double>>>&, vector<vector<complex<double>>>&, vector<vector<complex<double>>>&, vector<complex<double>>&, vector<complex<double>>&, map<string, int>&, double) {
        complex<double> conductance = 1.0 / (is_on ? Ron : Roff);
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
    void Circuit::build_ac_mna_matrix(vector<vector<complex<double>>>& A, vector<complex<double>>& z, double omega) {
        if (!components.empty() && !ground_node_exists) throw runtime_error("No ground node defined.");
        int N = idx_to_node_name.size();
        int M = voltage_source_list.size() + inductor_list.size() + vcvs_list.size() + ccvs_list.size();
        int system_size = N + M;
        if (system_size == 0) {
            A.clear(); z.clear(); return;
        }
        A.assign(system_size, vector<complex<double>>(system_size, 0.0));
        z.assign(system_size, 0.0);
        vector<vector<complex<double>>> G(N, vector<complex<double>>(N, 0.0)), B(N, vector<complex<double>>(M, 0.0)), C(M, vector<complex<double>>(N, 0.0)), D(M, vector<complex<double>>(M, 0.0));
        vector<complex<double>> J(N, 0.0), E(M, 0.0);
        map<string, int> m_map;
        int m_counter = 0;
        for (const auto& vs : voltage_source_list) m_map[vs->name] = m_counter++;
        for (const auto& l : inductor_list) m_map[l->name] = m_counter++;
        for (const auto& comp : components) comp->ac_stamp(*this, G, B, C, D, J, E, m_map, omega);
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
        tran_solved = false; dc_sweep_solved = false; ac_sweep_solved = false; phase_sweep_solved = false;
        tran_results.clear();
        if (t_step <= 0 || t_stop <= 0 || t_step > t_stop) throw runtime_error("Invalid transient parameters.");
        this->tran_t_stop = t_stop;
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
        tran_solved = false; dc_sweep_solved = false; ac_sweep_solved = false; phase_sweep_solved = false;
        dc_sweep_results.clear();
        if (inc == 0 || (end > start && inc < 0) || (end < start && inc > 0)) throw runtime_error("Invalid sweep parameters.");
        this->dc_start = start;
        this->dc_end = end;
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
    void Circuit::perform_ac_sweep_analysis(const string& src_name, double f_start, double f_end, int points_per_decade) {
        tran_solved = false; dc_sweep_solved = false; ac_sweep_solved = false; phase_sweep_solved = false;
        ac_sweep_results.clear();
        if (f_start <= 0 || f_end <= f_start || points_per_decade <= 0) throw runtime_error("Invalid AC sweep parameters.");
        int N = prepare_for_analysis();
        VoltageSource* sweep_src = nullptr;
        for (auto* vs : voltage_source_list) {
            if (vs->name == src_name && (vs->sourceType == VoltageSource::AC || vs->sourceType == VoltageSource::SINUSOIDAL)) {
                sweep_src = vs;
                break;
            }
        }
        if (!sweep_src) throw runtime_error("AC sweep source '" + src_name + "' not found or is not an AC/SIN source.");
        map<VoltageSource*, pair<double, double>> original_states;
        for (auto* vs : voltage_source_list) {
            if (vs->sourceType == VoltageSource::AC || vs->sourceType == VoltageSource::SINUSOIDAL) {
                original_states[vs] = {vs->amplitude, vs->phase_degrees};
                vs->amplitude = 0.0;
            }
        }
        sweep_src->amplitude = 1.0;
        sweep_src->phase_degrees = 0.0;
        int M = voltage_source_list.size() + inductor_list.size() + vcvs_list.size() + ccvs_list.size();
        if (N + M == 0) { ac_sweep_solved = true; return; }
        double freq_multiplier = pow(10.0, 1.0 / points_per_decade);
        for (double freq = f_start; freq <= f_end * (1 + 1e-9); freq *= freq_multiplier) {
            double omega = 2 * PI * freq;
            vector<vector<complex<double>>> A;
            vector<complex<double>> z;
            build_ac_mna_matrix(A, z, omega);
            vector<complex<double>> solution = complex_gaussian_elimination(A, z);
            ResultPoint result_at_freq;
            result_at_freq["frequency"] = freq;
            for (int i = 0; i < N; ++i) {
                string node_name = idx_to_node_name[i];
                complex<double> v = solution[i];
                double v_abs = abs(v);
                double v_mag_db = (v_abs > 1e-12) ? 20 * log10(v_abs) : -240.0;
                double v_phase = arg(v) * 180.0 / PI;
                result_at_freq["V(" + node_name + ")"] = v_mag_db;
                result_at_freq["V_mag(" + node_name + ")"] = v_abs;
                result_at_freq["V_mag_db(" + node_name + ")"] = v_mag_db;
                result_at_freq["V_phase(" + node_name + ")"] = v_phase;
            }
            map<string, int> m_map;
            int m_counter = 0;
            for (const auto& vs : voltage_source_list) m_map[vs->name] = m_counter++;
            for (const auto& l : inductor_list) m_map[l->name] = m_counter++;
            for (const auto& p : m_map) {
                string comp_name = p.first;
                complex<double> i = solution[N + p.second];
                double i_abs = abs(i);
                double i_mag_db = (i_abs > 1e-12) ? 20 * log10(i_abs) : -240.0;
                double i_phase = arg(i) * 180.0 / PI;
                result_at_freq["I(" + comp_name + ")"] = i_mag_db;
                result_at_freq["I_mag(" + comp_name + ")"] = i_abs;
                result_at_freq["I_mag_db(" + comp_name + ")"] = i_mag_db;
                result_at_freq["I_phase(" + comp_name + ")"] = i_phase;
            }
            ac_sweep_results.push_back(result_at_freq);
        }
        for (map<VoltageSource*, pair<double, double>>::const_iterator it = original_states.begin(); it != original_states.end(); ++it) {
            it->first->amplitude = it->second.first;
            it->first->phase_degrees = it->second.second;
        }
        ac_sweep_solved = true;
    }
    void Circuit::perform_phase_sweep_analysis(double freq, const string& src, double p_start, double p_end, double p_inc) {
        tran_solved = false; dc_sweep_solved = false; ac_sweep_solved = false; phase_sweep_solved = false;
        phase_sweep_results.clear();
        if (freq <= 0 || p_inc == 0) throw runtime_error("Invalid phase sweep parameters.");
        int N = prepare_for_analysis();
        VoltageSource* sweep_src = nullptr;
        for (auto* vs : voltage_source_list) {
            if (vs->name == src && (vs->sourceType == VoltageSource::AC || vs->sourceType == VoltageSource::SINUSOIDAL)) {
                sweep_src = vs;
                break;
            }
        }
        if (!sweep_src) throw runtime_error("AC sweep source not found or is not an AC/SIN source.");
        double orig_phase = sweep_src->phase_degrees;
        int M = voltage_source_list.size() + inductor_list.size();
        if (N + M == 0) { phase_sweep_solved = true; return; }
        double omega = 2 * PI * freq;
        for (double phase = p_start; (p_inc > 0) ? (phase <= p_end + abs(p_inc)/2.0) : (phase >= p_end - abs(p_inc)/2.0); phase += p_inc) {
            sweep_src->phase_degrees = phase;
            vector<vector<complex<double>>> A;
            vector<complex<double>> z;
            build_ac_mna_matrix(A, z, omega);
            vector<complex<double>> solution = complex_gaussian_elimination(A, z);
            ResultPoint result_at_phase;
            result_at_phase["phase"] = phase;
            for (int i = 0; i < N; ++i) {
                string node_name = idx_to_node_name[i];
                complex<double> v = solution[i];
                double v_mag = abs(v);
                double v_phase = arg(v) * 180.0 / PI;
                result_at_phase["V(" + node_name + ")"] = v_mag;
                result_at_phase["V_mag(" + node_name + ")"] = v_mag;
                result_at_phase["V_phase(" + node_name + ")"] = v_phase;
            }
            map<string, int> m_map;
            int m_counter = 0;
            for (const auto& vs : voltage_source_list) m_map[vs->name] = m_counter++;
            for (const auto& l : inductor_list) m_map[l->name] = m_counter++;
            for (const auto& p : m_map) {
                string comp_name = p.first;
                complex<double> i = solution[N + p.second];
                double i_mag = abs(i);
                double i_phase = arg(i) * 180.0 / PI;
                result_at_phase["I(" + comp_name + ")"] = i_mag;
                result_at_phase["I_mag(" + comp_name + ")"] = i_mag;
                result_at_phase["I_phase(" + comp_name + ")"] = i_phase;
            }
            phase_sweep_results.push_back(result_at_phase);
        }
        sweep_src->phase_degrees = orig_phase;
        phase_sweep_solved = true;
    }
}
unique_ptr<SpiceEngine::Circuit> last_simulated_circuit;
bool show_trace_dialog = false;
string new_trace_buffer = "";
int active_trace_dialog_field = -1;
vector<SDL_Color> plot_colors = {
        {255, 80, 80, 255}, {80, 255, 80, 255}, {80, 80, 255, 255},
        {255, 255, 80, 255}, {80, 255, 255, 255}, {255, 80, 255, 255}
};
bool isFileMenuOpen = false;
bool isComponentsMenuOpen = false;
bool isSimulateMenuOpen = false;
enum class InteractionMode {
    NONE, WIRING, PLACE_COMPONENT, DRAGGING_COMPONENT, DELETE_ITEM,
    EDITING_COMPONENT_VALUE, PLACE_LABEL, EDITING_LABEL_TEXT, PLACE_GND_LABEL, DIALOG_ACTIVE,
    PROBE_CURRENT, PROBE_VOLTAGE, PROBE_POWER
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
enum class DialogType { NONE, SAVE_AS, OPEN, TRANSIENT_ANALYSIS, DC_SWEEP_ANALYSIS, AC_SWEEP_ANALYSIS, PHASE_SWEEP_ANALYSIS };
DialogType activeDialogType = DialogType::NONE;
struct DialogField {
    string label, buffer;
    SDL_Rect input_rect;
};
vector<DialogField> dialogFields;
string dialogTitle;
int activeDialogFieldIndex = -1;
vector<string> dialog_file_list;
int dialog_file_list_scroll_offset = 0;
int dialog_hovered_file_index = -1;
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
void save_schematic(const string& filename_with_path) {
    ofstream outFile(filename_with_path);
    if (!outFile.is_open()) {
        cerr << "Error: Could not open file for saving: " << filename_with_path << endl;
        return;
    }
    for (const auto& comp : components) outFile << "COMPONENT," << component_type_to_string(comp.type) << "," << comp.id << "," << comp.node1.x << "," << comp.node1.y << "," << comp.node2.x << "," << comp.node2.y << "," << comp.rotation_angle << "," << comp.value << endl;
    for (const auto& wire : wires) outFile << "WIRE," << wire.start.x << "," << wire.start.y << "," << wire.end.x << "," << wire.end.y << endl;
    for (const auto& label : labels) outFile << "LABEL," << label.text << "," << label.position.x << "," << label.position.y << "," << label.rotation_angle << endl;
    outFile.close();
    size_t last_slash = filename_with_path.find_last_of("/\\");
    currentSchematicFileName = (last_slash == string::npos) ? filename_with_path : filename_with_path.substr(last_slash + 1);
    schematicModified = false;
    cout << "Schematic saved to " << filename_with_path << endl;
}
void load_schematic(const string& filename_with_path) {
    ifstream inFile(filename_with_path);
    if (!inFile.is_open()) {
        cerr << "Error: Could not open file for loading: " << filename_with_path << endl;
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
    size_t last_slash = filename_with_path.find_last_of("/\\");
    currentSchematicFileName = (last_slash == string::npos) ? filename_with_path : filename_with_path.substr(last_slash + 1);
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
map<string, string> resolve_node_names(string& out_ground_node_name) {
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
    out_ground_node_name = "";
    for (const auto& label : labels) {
        if (label.text == "GND") {
            string key = point_to_key(label.position);
            if (point_to_set_id.count(key)) {
                int root = find_set(point_to_set_id.at(key));
                out_ground_node_name = set_to_node_name[root];
                break;
            }
        }
    }
    map<string, string> point_to_node_name_map;
    for (const auto& pair : point_to_set_id) {
        string point_key = pair.first;
        int root = find_set(pair.second);
        point_to_node_name_map[point_key] = set_to_node_name.at(root);
    }
    return point_to_node_name_map;
}
void run_simulation_in_terminal(DialogType analysis_type, const vector<DialogField>& fields, unique_ptr<SpiceEngine::Circuit>& circuit_out) {
    cout << "\n--- Preparing Simulation ---" << endl;
    circuit_out = make_unique<SpiceEngine::Circuit>();
    component_id_to_nodes.clear();
    string ground_node_name;
    auto point_to_node_name_map = resolve_node_names(ground_node_name);
    if (!ground_node_name.empty()) {
        circuit_out->set_ground_node(ground_node_name);
    }
    auto point_to_key = [](SDL_Point p) { return to_string(p.x) + "," + to_string(p.y); };
    for (const auto& comp : components) {
        string n1 = point_to_node_name_map.at(point_to_key(comp.node1));
        string n2 = point_to_node_name_map.at(point_to_key(comp.node2));
        component_id_to_nodes[comp.id] = {n1, n2};
        switch (comp.type) {
            case ComponentType::RESISTOR: circuit_out->add_component(make_unique<SpiceEngine::Resistor>(comp.id, n1, n2, comp.value)); break;
            case ComponentType::CAPACITOR: circuit_out->add_component(make_unique<SpiceEngine::Capacitor>(comp.id, n1, n2, comp.value)); break;
            case ComponentType::INDUCTOR: circuit_out->add_component(make_unique<SpiceEngine::Inductor>(comp.id, n1, n2, comp.value)); break;
            case ComponentType::DIODE: circuit_out->add_component(make_unique<SpiceEngine::Diode>(comp.id, n1, n2, comp.value)); break;
            case ComponentType::DC_CURRENT_SOURCE: circuit_out->add_component(make_unique<SpiceEngine::CurrentSource>(comp.id, n1, n2, comp.value)); break;
            case ComponentType::DC_VOLTAGE_SOURCE:
            case ComponentType::AC_VOLTAGE_SOURCE: {
                vector<string> params;
                string temp_val = comp.value;
                if (to_lower_util(temp_val).find("sin") != string::npos || to_lower_util(temp_val).find("pulse") != string::npos) {
                    params.push_back(temp_val);
                } else {
                    stringstream ss(temp_val);
                    string token;
                    while(ss >> token) {
                        params.push_back(token);
                    }
                }
                circuit_out->add_component(make_unique<SpiceEngine::VoltageSource>(comp.id, n1, n2, params));
                break;
            }
            default: break;
        }
    }
    try {
        plotted_variables.clear();
        vector<string> requested_vars;
        if (analysis_type == DialogType::TRANSIENT_ANALYSIS) {
            double t_step = parse_value_with_metric_prefix_util(fields[0].buffer);
            double t_stop = parse_value_with_metric_prefix_util(fields[1].buffer);
            stringstream ss(fields[2].buffer);
            string var;
            while(ss >> var) { requested_vars.push_back(var); }
            circuit_out->perform_transient_analysis(t_step, t_stop);
            cout << left << setw(15) << "time (s)";
            for (const auto& v : requested_vars) cout << setw(15) << v;
            cout << endl;
            for (const auto& result_point : circuit_out->tran_results) {
                cout << fixed << setprecision(6) << setw(15) << result_point.at("time");
                for (const auto& v : requested_vars) {
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
            while(ss >> var) { requested_vars.push_back(var); }
            circuit_out->perform_dc_sweep_analysis(src_name, start, end, inc);
            string header_src_name = circuit_out->dc_sweep_source_name;
            char src_name_char = to_lower_util(header_src_name)[0];
            if (src_name_char == 'v') header_src_name += " (V)";
            else if (src_name_char == 'i') header_src_name += " (A)";
            cout << left << setw(15) << header_src_name;
            for (const auto& v : requested_vars) cout << setw(15) << v;
            cout << endl;
            for (const auto& result_point : circuit_out->dc_sweep_results) {
                cout << fixed << setprecision(6) << setw(15) << result_point.at(circuit_out->dc_sweep_source_name);
                for (const auto& v : requested_vars) {
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
        } else if (analysis_type == DialogType::AC_SWEEP_ANALYSIS) {
            string src_name = fields[0].buffer;
            double f_start = parse_value_with_metric_prefix_util(fields[1].buffer);
            double f_end = parse_value_with_metric_prefix_util(fields[2].buffer);
            int points = stoi(fields[3].buffer);
            stringstream ss(fields[4].buffer);
            string var;
            while(ss >> var) { requested_vars.push_back(var); }
            circuit_out->perform_ac_sweep_analysis(src_name, f_start, f_end, points);
            cout << left << setw(15) << "frequency (Hz)";
            for (const auto& v : requested_vars) cout << setw(20) << v;
            cout << endl;
            for (const auto& result_point : circuit_out->ac_sweep_results) {
                cout << scientific << setprecision(6) << setw(15) << result_point.at("frequency");
                for (const auto& v : requested_vars) {
                    if (result_point.count(v)) {
                        cout << fixed << setprecision(6) << setw(20) << result_point.at(v);
                    } else {
                        cout << setw(20) << "N/A";
                    }
                }
                cout << endl;
            }
        } else if (analysis_type == DialogType::PHASE_SWEEP_ANALYSIS) {
            double freq = parse_value_with_metric_prefix_util(fields[0].buffer);
            string src = fields[1].buffer;
            double p_start = parse_value_with_metric_prefix_util(fields[2].buffer);
            double p_end = parse_value_with_metric_prefix_util(fields[3].buffer);
            double p_inc = parse_value_with_metric_prefix_util(fields[4].buffer);
            stringstream ss(fields[5].buffer);
            string var;
            while(ss >> var) { requested_vars.push_back(var); }
            circuit_out->perform_phase_sweep_analysis(freq, src, p_start, p_end, p_inc);
            cout << left << setw(15) << "phase (deg)";
            for (const auto& v : requested_vars) cout << setw(20) << v;
            cout << endl;
            for (const auto& result_point : circuit_out->phase_sweep_results) {
                cout << fixed << setprecision(6) << setw(15) << result_point.at("phase");
                for (const auto& v : requested_vars) {
                    if (result_point.count(v)) {
                        cout << setw(20) << result_point.at(v);
                    } else {
                        cout << setw(20) << "N/A";
                    }
                }
                cout << endl;
            }
        }
        for(size_t i = 0; i < requested_vars.size(); ++i) {
            plotted_variables.push_back({requested_vars[i], plot_colors[i % plot_colors.size()]});
        }
    } catch (const exception& e) {
        cerr << "SIMULATION ERROR: " << e.what() << endl;
        circuit_out.reset();
    }
    cout << "--- Simulation Finished ---" << endl;
}
void close_dialog() {
    activeDialogType = DialogType::NONE;
    dialogFields.clear();
    dialog_file_list.clear();
    dialogTitle = "";
    activeDialogFieldIndex = -1;
    currentInteractionMode = InteractionMode::NONE;
    SDL_StopTextInput();
}
void handle_dialog_ok() {
    if (activeDialogType == DialogType::SAVE_AS) {
        if (!dialogFields.empty() && !dialogFields[0].buffer.empty()) {
            string filename = dialogFields[0].buffer;
            if (filename.length() < 4 || to_lower_util(filename.substr(filename.length() - 4)) != ".txt") {
                filename += ".txt";
            }
            save_schematic(SCHEMATICS_PATH + filename);
        }
    }
    else if (activeDialogType == DialogType::TRANSIENT_ANALYSIS || activeDialogType == DialogType::DC_SWEEP_ANALYSIS || activeDialogType == DialogType::AC_SWEEP_ANALYSIS || activeDialogType == DialogType::PHASE_SWEEP_ANALYSIS) {
        run_simulation_in_terminal(activeDialogType, dialogFields, last_simulated_circuit);
        if (last_simulated_circuit && (last_simulated_circuit->tran_solved || last_simulated_circuit->dc_sweep_solved || last_simulated_circuit->ac_sweep_solved || last_simulated_circuit->phase_sweep_solved)) {
            currentAppState = AppState::RESULTS_VIEW;
            cursor_index = 0;
        }
    }
    close_dialog();
}
void setup_dialog(DialogType type, const string& title, const vector<string>& labels) {
    activeDialogType = type;
    dialogTitle = title;
    dialogFields.clear();
    dialog_file_list.clear();
    if (type == DialogType::OPEN) {
        dialog_file_list = get_schematic_files();
        dialog_file_list_scroll_offset = 0;
    } else {
        for (const auto& label : labels) dialogFields.push_back({ label, "" });
        activeDialogFieldIndex = 0;
        SDL_StartTextInput();
    }
    currentInteractionMode = InteractionMode::DIALOG_ACTIVE;
}
void render_dialog(SDL_Renderer* renderer, TTF_Font* font) {
    if (activeDialogType == DialogType::OPEN) {
        const int dialog_w = 450, dialog_h = 400;
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
        const int item_height = 25;
        const int list_height = dialog_h - 120;
        const int max_visible_items = list_height / item_height;
        SDL_Rect list_view_rect = {panelRect.x + 10, current_y, panelRect.w - 20, list_height};
        SDL_RenderSetClipRect(renderer, &list_view_rect);
        for (size_t i = 0; i < dialog_file_list.size(); ++i) {
            int displayed_index = i - dialog_file_list_scroll_offset;
            if (displayed_index < 0 || displayed_index >= max_visible_items) continue;
            SDL_Rect item_rect = { panelRect.x + 20, current_y + displayed_index * item_height, panelRect.w - 40, item_height };
            if (dialog_hovered_file_index == (int)i) {
                SDL_SetRenderDrawColor(renderer, 0x77, 0x77, 0x99, 0xFF);
                SDL_RenderFillRect(renderer, &item_rect);
            }
            SDL_Surface* fileSurface = TTF_RenderText_Solid(font, dialog_file_list[i].c_str(), textColor);
            if (fileSurface) {
                SDL_Texture* fileTexture = SDL_CreateTextureFromSurface(renderer, fileSurface);
                SDL_Rect fileRect = { item_rect.x + 5, item_rect.y + (item_rect.h - fileSurface->h) / 2, fileSurface->w, fileSurface->h };
                SDL_RenderCopy(renderer, fileTexture, nullptr, &fileRect);
                SDL_FreeSurface(fileSurface);
                SDL_DestroyTexture(fileTexture);
            }
        }
        SDL_RenderSetClipRect(renderer, nullptr);
        Button cancelButton("Cancel", panelRect.x + panelRect.w / 2 - 40, panelRect.y + panelRect.h - 45, 80, 30, close_dialog);
        int mx, my;
        SDL_GetMouseState(&mx, &my);
        cancelButton.m_is_hovered = SDL_PointInRect(new const SDL_Point{ mx, my }, &cancelButton.m_position);
        cancelButton.render(renderer, font);
    } else {
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
            SDL_SetRenderDrawColor(renderer, (activeDialogFieldIndex == (int)i) ? 0x00 : 0x33, (activeDialogFieldIndex == (int)i) ? 0x00 : 0x33, (activeDialogFieldIndex == (int)i) ? 0x44 : 0x33, 0xFF);
            SDL_RenderFillRect(renderer, &inputRect);
            SDL_SetRenderDrawColor(renderer, 0xAA, 0xAA, 0xAA, 0xFF);
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
}
SDL_Point map_point_to_screen(double data_x, double data_y, const SDL_Rect& graph_area, double x_min, double x_max, double y_min, double y_max, bool log_x) {
    if (x_max == x_min || y_max == y_min) {
        return { graph_area.x, graph_area.y };
    }
    int screen_x;
    if (log_x) {
        if(data_x <= 0) data_x = x_min;
        screen_x = graph_area.x + static_cast<int>(((log10(data_x) - log10(x_min)) / (log10(x_max) - log10(x_min))) * graph_area.w);
    } else {
        screen_x = graph_area.x + static_cast<int>(((data_x - x_min) / (x_max - x_min)) * graph_area.w);
    }
    int screen_y = graph_area.y + graph_area.h - static_cast<int>(((data_y - y_min) / (y_max - y_min)) * graph_area.h);
    return { screen_x, screen_y };
}
void render_text(SDL_Renderer* renderer, TTF_Font* font, const string& text, int x, int y, SDL_Color color, bool center_x = false, bool center_y = false) {
    if (!font || text.empty()) return;
    SDL_Surface* surface = TTF_RenderText_Solid(font, text.c_str(), color);
    if (!surface) return;
    SDL_Texture* texture = SDL_CreateTextureFromSurface(renderer, surface);
    if (!texture) {
        SDL_FreeSurface(surface);
        return;
    }
    SDL_Rect dest_rect = { x, y, surface->w, surface->h };
    if (center_x) dest_rect.x -= surface->w / 2;
    if (center_y) dest_rect.y -= surface->h / 2;
    SDL_RenderCopy(renderer, texture, nullptr, &dest_rect);
    SDL_DestroyTexture(texture);
    SDL_FreeSurface(surface);
}
double evaluate_expression(const string& expr, const SpiceEngine::ResultPoint& data_point) {
    string trimmed_expr = trim_string_util(expr);
    if (trimmed_expr.length() > 3 && toupper(trimmed_expr[0]) == 'P' && trimmed_expr[1] == '(' && trimmed_expr.back() == ')') {
        string comp_id = trimmed_expr.substr(2, trimmed_expr.length() - 3);
        if (component_id_to_nodes.count(comp_id)) {
            pair<string, string> nodes = component_id_to_nodes.at(comp_id);
            string v1_var = "V(" + nodes.first + ")";
            string v2_var = "V(" + nodes.second + ")";
            string i_var = "I(" + comp_id + ")";
            double v1 = data_point.count(v1_var) ? data_point.at(v1_var) : 0.0;
            double v2 = data_point.count(v2_var) ? data_point.at(v2_var) : 0.0;
            double i = data_point.count(i_var) ? data_point.at(i_var) : NAN;
            if (isnan(i)) return NAN;
            return (v1 - v2) * i;
        }
        return NAN;
    }
    string current_token;
    vector<string> tokens;
    vector<char> ops;
    for (char c : trimmed_expr) {
        if (c == '+' || c == '-') {
            if (!current_token.empty()) {
                tokens.push_back(trim_string_util(current_token));
                current_token.clear();
            }
            ops.push_back(c);
        } else {
            current_token += c;
        }
    }
    if (!current_token.empty()) {
        tokens.push_back(trim_string_util(current_token));
    }
    if (tokens.empty()) return NAN;
    double result = 0.0;
    if (data_point.count(tokens[0])) {
        result = data_point.at(tokens[0]);
    } else {
        return NAN;
    }
    for (size_t i = 0; i < ops.size(); ++i) {
        if (i + 1 < tokens.size()) {
            if (data_point.count(tokens[i+1])) {
                double val = data_point.at(tokens[i+1]);
                if (ops[i] == '+') result += val;
                else if (ops[i] == '-') result -= val;
            } else {
                return NAN;
            }
        } else {
            return NAN;
        }
    }
    return result;
}
void render_trace_dialog(SDL_Renderer* renderer, TTF_Font* font);
void render_scale_dialog(SDL_Renderer* renderer, TTF_Font* font);
void render_results_view(SDL_Renderer* renderer, TTF_Font* font, vector<Button>& results_buttons, const vector<ComponentMenuItem>& componentMenuItems) {
    SDL_Rect controlPanel = { 0, 0, 220, SCREEN_HEIGHT };
    SDL_SetRenderDrawColor(renderer, 0x33, 0x33, 0x33, 0xFF);
    SDL_RenderFillRect(renderer, &controlPanel);
    const int margin_top = 50, margin_bottom = 50, margin_left = 70, margin_right = 30;
    SDL_Rect graphArea = { controlPanel.w + margin_left, margin_top, SCREEN_WIDTH - controlPanel.w - margin_left - margin_right, SCREEN_HEIGHT - margin_top - margin_bottom };
    SDL_SetRenderDrawColor(renderer, 0x11, 0x11, 0x11, 0xFF);
    SDL_RenderFillRect(renderer, &graphArea);
    render_text(renderer, font, "Simulation Results", controlPanel.w / 2, 20, { 255, 255, 255, 255 }, true);
    for (auto& button : results_buttons) {
        button.render(renderer, font);
    }
    if (!last_simulated_circuit || plotted_variables.empty()) {
        render_text(renderer, font, "No data to display", graphArea.x + graphArea.w / 2, graphArea.y + graphArea.h / 2, { 150, 150, 150, 255 }, true, true);
        SDL_SetRenderDrawColor(renderer, 0x88, 0x88, 0x88, 0xFF);
        SDL_RenderDrawRect(renderer, &graphArea);
        if(show_trace_dialog) render_trace_dialog(renderer, font);
        if(show_scale_dialog) render_scale_dialog(renderer, font);
        return;
    }
    const bool is_tran = last_simulated_circuit->tran_solved;
    const bool is_dc = last_simulated_circuit->dc_sweep_solved;
    const bool is_ac = last_simulated_circuit->ac_sweep_solved;
    const bool is_phase = last_simulated_circuit->phase_sweep_solved;
    const auto& results = is_tran ? last_simulated_circuit->tran_results :
                          is_dc ? last_simulated_circuit->dc_sweep_results :
                          is_ac ? last_simulated_circuit->ac_sweep_results :
                          last_simulated_circuit->phase_sweep_results;
    if (results.empty()) {
        render_text(renderer, font, "Simulation produced no data points", graphArea.x + graphArea.w / 2, graphArea.y + graphArea.h / 2, { 150, 150, 150, 255 }, true, true);
        SDL_SetRenderDrawColor(renderer, 0x88, 0x88, 0x88, 0xFF);
        SDL_RenderDrawRect(renderer, &graphArea);
        if(show_trace_dialog) render_trace_dialog(renderer, font);
        if(show_scale_dialog) render_scale_dialog(renderer, font);
        return;
    }
    bool is_log_x = is_ac;
    if (!results.empty() && results.front().count("frequency")) {
        is_log_x = true;
    }
    double y_min_auto = DBL_MAX, y_max_auto = -DBL_MAX;
    for (const auto& var : plotted_variables) {
        for (const auto& point : results) {
            double val = evaluate_expression(var.name, point);
            if (!isnan(val)) {
                if (val < y_min_auto) y_min_auto = val;
                if (val > y_max_auto) y_max_auto = val;
            }
        }
    }
    if (y_min_auto > y_max_auto) {
        render_text(renderer, font, "Variable not found in results", graphArea.x + graphArea.w / 2, graphArea.y + graphArea.h / 2, { 255, 100, 100, 255 }, true, true);
        SDL_SetRenderDrawColor(renderer, 0x88, 0x88, 0x88, 0xFF);
        SDL_RenderDrawRect(renderer, &graphArea);
        if(show_trace_dialog) render_trace_dialog(renderer, font);
        if(show_scale_dialog) render_scale_dialog(renderer, font);
        return;
    }
    double x_min, x_max;
    if (is_tran) { x_min = 0.0; x_max = last_simulated_circuit->tran_t_stop; }
    else if (is_dc) { x_min = last_simulated_circuit->dc_start; x_max = last_simulated_circuit->dc_end; }
    else if (is_ac) { x_min = results.front().at("frequency"); x_max = results.back().at("frequency"); }
    else { x_min = results.front().at("phase"); x_max = results.back().at("phase"); }
    if (!manual_scale.x_min_str.empty()) x_min = parse_value_with_metric_prefix_util(manual_scale.x_min_str);
    if (!manual_scale.x_max_str.empty()) x_max = parse_value_with_metric_prefix_util(manual_scale.x_max_str);
    double y_min = y_min_auto;
    double y_max = y_max_auto;
    double y_range = y_max - y_min;
    if (abs(y_range) < 1e-9) {
        y_range = (abs(y_min) > 1e-9) ? abs(y_min) : 1.0;
        y_min -= y_range * 0.5;
        y_max += y_range * 0.5;
    }
    y_min -= y_range * 0.1;
    y_max += y_range * 0.1;
    if (!manual_scale.y_min_str.empty()) y_min = parse_value_with_metric_prefix_util(manual_scale.y_min_str);
    if (!manual_scale.y_max_str.empty()) y_max = parse_value_with_metric_prefix_util(manual_scale.y_max_str);
    y_range = y_max - y_min;
    const int num_grid_lines = 10;
    SDL_SetRenderDrawColor(renderer, 0x40, 0x40, 0x40, 0xFF);
    for (int i = 1; i < num_grid_lines; ++i) {
        double x_val;
        if(is_log_x) {
            x_val = x_min * pow(x_max/x_min, (double)i/num_grid_lines);
        } else {
            x_val = x_min + (x_max - x_min) * i / num_grid_lines;
        }
        int x_screen = map_point_to_screen(x_val, 0, graphArea, x_min, x_max, y_min, y_max, is_log_x).x;
        SDL_RenderDrawLine(renderer, x_screen, graphArea.y, x_screen, graphArea.y + graphArea.h);
        render_text(renderer, font, format_value_with_metric_prefix(x_val), x_screen, graphArea.y + graphArea.h + 5, {200, 200, 200, 255}, true);
    }
    for (int i = 1; i < num_grid_lines; ++i) {
        double y_val = y_min + y_range * i / num_grid_lines;
        int y_screen = graphArea.y + graphArea.h - (graphArea.h * i / num_grid_lines);
        SDL_RenderDrawLine(renderer, graphArea.x, y_screen, graphArea.x + graphArea.w, y_screen);
        render_text(renderer, font, format_value_with_metric_prefix(y_val), graphArea.x - 10, y_screen, {200, 200, 200, 255}, true, true);
    }
    SDL_SetRenderDrawColor(renderer, 0x88, 0x88, 0x88, 0xFF);
    SDL_RenderDrawRect(renderer, &graphArea);
    render_text(renderer, font, format_value_with_metric_prefix(y_max), graphArea.x - 10, graphArea.y, { 200, 200, 200, 255 }, true, true);
    render_text(renderer, font, format_value_with_metric_prefix(y_min), graphArea.x - 10, graphArea.y + graphArea.h, { 200, 200, 200, 255 }, true, true);
    render_text(renderer, font, format_value_with_metric_prefix(x_min), graphArea.x, graphArea.y + graphArea.h + 5, { 200, 200, 200, 255 }, true);
    render_text(renderer, font, format_value_with_metric_prefix(x_max), graphArea.x + graphArea.w, graphArea.y + graphArea.h + 5, { 200, 200, 200, 255 }, true);
    string x_axis_key = is_tran ? "time" : is_dc ? last_simulated_circuit->dc_sweep_source_name : is_ac ? "frequency" : "phase";
    int legend_y = 270;
    for (size_t i = 0; i < plotted_variables.size(); ++i) {
        const auto& var = plotted_variables[i];
        SDL_Rect legend_color_box = { 10, legend_y, 15, 15 };
        SDL_SetRenderDrawColor(renderer, var.color.r, var.color.g, var.color.b, 255);
        SDL_RenderFillRect(renderer, &legend_color_box);
        render_text(renderer, font, var.name, 35, legend_y, { 255, 255, 255, 255 });
        legend_y += 25;
        SDL_Point prev_point = { -1, -1 };
        for (const auto& point_data : results) {
            if (point_data.count(x_axis_key)) {
                double x_val = point_data.at(x_axis_key);
                double y_val = evaluate_expression(var.name, point_data);
                if (!isnan(y_val)) {
                    SDL_Point current_point = map_point_to_screen(x_val, y_val, graphArea, x_min, x_max, y_min, y_max, is_log_x);
                    if (prev_point.x != -1) {
                        thickLineRGBA(renderer, prev_point.x, prev_point.y, current_point.x, current_point.y, 2, var.color.r, var.color.g, var.color.b, 255);
                    }
                    prev_point = current_point;
                } else {
                    prev_point = {-1, -1};
                }
            }
        }
    }
    if (cursor_index >= 0 && cursor_index < (int)results.size()) {
        const auto& cursor_data_point = results[cursor_index];
        double x_cursor_val = cursor_data_point.at(x_axis_key);
        SDL_Point cursor_screen_pos = map_point_to_screen(x_cursor_val, 0, graphArea, x_min, x_max, y_min, y_max, is_log_x);
        SDL_SetRenderDrawColor(renderer, 255, 255, 0, 150);
        SDL_RenderDrawLine(renderer, cursor_screen_pos.x, graphArea.y, cursor_screen_pos.x, graphArea.y + graphArea.h);
        int info_y = legend_y + 20;
        render_text(renderer, font, "Cursor Data:", 110, info_y, {255, 255, 0, 255}, true);
        info_y += 25;
        render_text(renderer, font, x_axis_key + ": " + format_value_with_metric_prefix(x_cursor_val), 10, info_y, {255, 255, 255, 255});
        info_y += 20;
        for (const auto& var : plotted_variables) {
            double y_cursor_val = evaluate_expression(var.name, cursor_data_point);
            if (!isnan(y_cursor_val)) {
                SDL_Point trace_cursor_pos = map_point_to_screen(x_cursor_val, y_cursor_val, graphArea, x_min, x_max, y_min, y_max, is_log_x);
                filledCircleRGBA(renderer, trace_cursor_pos.x, trace_cursor_pos.y, 5, var.color.r, var.color.g, var.color.b, 255);
                string info_text = var.name + ": " + format_value_with_metric_prefix(y_cursor_val);
                render_text(renderer, font, info_text, 10, info_y, var.color);
                info_y += 20;
                if (cursor_index > 0 && cursor_index < (int)results.size() - 1) {
                    const auto& prev_point = results[cursor_index - 1];
                    const auto& next_point = results[cursor_index + 1];
                    double y_prev = evaluate_expression(var.name, prev_point);
                    double y_next = evaluate_expression(var.name, next_point);
                    double x_prev = prev_point.at(x_axis_key);
                    double x_next = next_point.at(x_axis_key);
                    if (!isnan(y_prev) && !isnan(y_next) && (x_next - x_prev) != 0) {
                        double slope = (y_next - y_prev) / (x_next - x_prev);
                        string slope_units = string(" / ") + (is_ac ? "Hz" : is_dc ? "V" : is_phase ? "deg" : "s");
                        string slope_text = "Slope: " + format_value_with_metric_prefix(slope) + slope_units;
                        render_text(renderer, font, slope_text, 10, info_y, var.color);
                        info_y += 20;
                    }
                }
            }
        }
    }
    if (currentInteractionMode == InteractionMode::PROBE_CURRENT || currentInteractionMode == InteractionMode::PROBE_VOLTAGE || currentInteractionMode == InteractionMode::PROBE_POWER) {
        SDL_SetRenderDrawColor(renderer, 0, 0, 0, 150);
        SDL_RenderFillRect(renderer, nullptr);
        draw_schematic_elements(renderer, font, componentMenuItems);
        int mx, my;
        SDL_GetMouseState(&mx, &my);
        if (currentInteractionMode == InteractionMode::PROBE_CURRENT) {
            render_text(renderer, font, "Click on a component to probe its current. Press ESC to cancel.", SCREEN_WIDTH / 2, 20, {255, 255, 0, 255}, true);
            circleRGBA(renderer, mx, my, 15, 255, 0, 0, 255);
            lineRGBA(renderer, mx - 20, my, mx + 20, my, 255, 0, 0, 255);
            lineRGBA(renderer, mx, my - 20, mx, my + 20, 255, 0, 0, 255);
        } else if (currentInteractionMode == InteractionMode::PROBE_VOLTAGE) {
            render_text(renderer, font, "Click on a node to probe its voltage. Press ESC to cancel.", SCREEN_WIDTH / 2, 20, {255, 255, 0, 255}, true);
            thickLineRGBA(renderer, mx, my, mx, my - 20, 3, 0, 255, 0, 255);
            thickLineRGBA(renderer, mx, my, mx - 14, my + 14, 3, 0, 255, 0, 255);
            thickLineRGBA(renderer, mx, my, mx + 14, my + 14, 3, 0, 255, 0, 255);
        } else {
            render_text(renderer, font, "Click on a component to probe its power. Press ESC to cancel.", SCREEN_WIDTH / 2, 20, {255, 255, 0, 255}, true);
            Sint16 vx[] = { (Sint16)(mx), (Sint16)(mx-8), (Sint16)(mx+2), (Sint16)(mx-2), (Sint16)(mx+8), (Sint16)(mx-2), (Sint16)(mx) };
            Sint16 vy[] = { (Sint16)(my-15), (Sint16)(my-5), (Sint16)(my-5), (Sint16)(my+5), (Sint16)(my+5), (Sint16)(my+15), (Sint16)(my+5) };
            filledPolygonRGBA(renderer, vx, vy, 7, 255, 215, 0, 255);
        }
    }
    if(show_trace_dialog) render_trace_dialog(renderer, font);
    if(show_scale_dialog) render_scale_dialog(renderer, font);
}
void render_trace_dialog(SDL_Renderer* renderer, TTF_Font* font) {
    const int dialog_w = 400;
    const int dialog_h = 150 + plotted_variables.size() * 40;
    SDL_Rect panelRect = { SCREEN_WIDTH / 2 - dialog_w / 2, SCREEN_HEIGHT / 2 - dialog_h / 2, dialog_w, dialog_h };
    SDL_SetRenderDrawColor(renderer, 0x55, 0x55, 0x55, 0xFF);
    SDL_RenderFillRect(renderer, &panelRect);
    SDL_SetRenderDrawColor(renderer, 0x88, 0x88, 0x88, 0xFF);
    SDL_RenderDrawRect(renderer, &panelRect);
    render_text(renderer, font, "Manage Traces", panelRect.x + dialog_w / 2, panelRect.y + 20, {255, 255, 255, 255}, true);
    int current_y = panelRect.y + 50;
    for(size_t i = 0; i < plotted_variables.size(); ++i) {
        SDL_Rect color_box = {panelRect.x + 20, current_y, 30, 30};
        SDL_SetRenderDrawColor(renderer, plotted_variables[i].color.r, plotted_variables[i].color.g, plotted_variables[i].color.b, 255);
        SDL_RenderFillRect(renderer, &color_box);
        render_text(renderer, font, plotted_variables[i].name, panelRect.x + 60, current_y + 15, {255, 255, 255, 255}, false, true);
        Button delete_button("Delete", panelRect.x + dialog_w - 80, current_y, 60, 30, [i](){
            if (i < plotted_variables.size()) {
                plotted_variables.erase(plotted_variables.begin() + i);
            }
        });
        delete_button.render(renderer, font);
        current_y += 40;
    }
    render_text(renderer, font, "Add:", panelRect.x + 20, current_y + 15, {255, 255, 255, 255}, false, true);
    SDL_Rect input_rect = {panelRect.x + 60, current_y, dialog_w - 150, 30};
    SDL_SetRenderDrawColor(renderer, (active_trace_dialog_field == 0) ? 0x00 : 0x33, (active_trace_dialog_field == 0) ? 0x00 : 0x33, (active_trace_dialog_field == 0) ? 0x44 : 0x33, 0xFF);
    SDL_RenderFillRect(renderer, &input_rect);
    render_text(renderer, font, new_trace_buffer, input_rect.x + 5, input_rect.y + 15, {255, 255, 255, 255}, false, true);
    Button add_button("Add", input_rect.x + input_rect.w + 10, current_y, 60, 30, [](){
        if(!new_trace_buffer.empty()) {
            plotted_variables.push_back({new_trace_buffer, plot_colors[plotted_variables.size() % plot_colors.size()]});
            new_trace_buffer = "";
        }
    });
    add_button.render(renderer, font);
    Button close_button("Close", panelRect.x + dialog_w/2 - 40, panelRect.y + dialog_h - 45, 80, 30, [](){ show_trace_dialog = false; SDL_StopTextInput(); });
    close_button.render(renderer, font);
}
void render_scale_dialog(SDL_Renderer* renderer, TTF_Font* font) {
    const int dialog_w = 350;
    const int dialog_h = 280;
    SDL_Rect panelRect = { SCREEN_WIDTH / 2 - dialog_w / 2, SCREEN_HEIGHT / 2 - dialog_h / 2, dialog_w, dialog_h };
    SDL_SetRenderDrawColor(renderer, 0x55, 0x55, 0x55, 0xFF);
    SDL_RenderFillRect(renderer, &panelRect);
    SDL_SetRenderDrawColor(renderer, 0x88, 0x88, 0x88, 0xFF);
    SDL_RenderDrawRect(renderer, &panelRect);
    render_text(renderer, font, "Adjust Scale", panelRect.x + dialog_w / 2, panelRect.y + 20, {255, 255, 255, 255}, true);
    vector<tuple<string, string*, int>> fields = {
            make_tuple("X-Min:", &manual_scale.x_min_str, 0),
            make_tuple("X-Max:", &manual_scale.x_max_str, 1),
            make_tuple("Y-Min:", &manual_scale.y_min_str, 2),
            make_tuple("Y-Max:", &manual_scale.y_max_str, 3)
    };
    int current_y = panelRect.y + 50;
    for (const auto& field : fields) {
        render_text(renderer, font, get<0>(field), panelRect.x + 20, current_y + 15, {255, 255, 255, 255}, false, true);
        SDL_Rect input_rect = {panelRect.x + 80, current_y, dialog_w - 100, 30};
        SDL_SetRenderDrawColor(renderer, (active_trace_dialog_field == get<2>(field)) ? 0x00 : 0x33, (active_trace_dialog_field == get<2>(field)) ? 0x00 : 0x33, (active_trace_dialog_field == get<2>(field)) ? 0x44 : 0x33, 0xFF);
        SDL_RenderFillRect(renderer, &input_rect);
        render_text(renderer, font, *get<1>(field), input_rect.x + 5, input_rect.y + 15, {255, 255, 255, 255}, false, true);
        current_y += 40;
    }
    Button close_button("Apply", panelRect.x + dialog_w/2 - 40, panelRect.y + dialog_h - 45, 80, 30, [](){ show_scale_dialog = false; SDL_StopTextInput(); });
    close_button.render(renderer, font);
}
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
    fileMenuButtons.emplace_back("Open...", 10, 70, 180, 30, []() { setup_dialog(DialogType::OPEN, "Open Schematic", {}); isFileMenuOpen = false; });
    fileMenuButtons.emplace_back("Save", 10, 100, 180, 30, []() { if (currentSchematicFileName.empty()) { setup_dialog(DialogType::SAVE_AS, "Save Schematic As", { "Filename" }); } else { save_schematic(SCHEMATICS_PATH + currentSchematicFileName); } isFileMenuOpen = false; });
    fileMenuButtons.emplace_back("Save As...", 10, 130, 180, 30, []() { setup_dialog(DialogType::SAVE_AS, "Save Schematic As", { "Filename" }); isFileMenuOpen = false; });
    bool quit_flag = false;
    fileMenuButtons.emplace_back("Exit", 10, 160, 180, 30, [&quit_flag]() { quit_flag = true; });
    vector<Button> simulateMenuButtons;
    simulateMenuButtons.emplace_back("DC Sweep Analysis", 100, 40, 180, 30, []() { setup_dialog(DialogType::DC_SWEEP_ANALYSIS, "DC Sweep Analysis", { "Source Name", "Start Value", "End Value", "Increment", "Wanted Values" }); isSimulateMenuOpen = false; });
    simulateMenuButtons.emplace_back("Transient Analysis", 100, 70, 180, 30, []() { setup_dialog(DialogType::TRANSIENT_ANALYSIS, "Transient Analysis", { "Tstep", "Tstop", "Wanted Values" }); isSimulateMenuOpen = false; });
    simulateMenuButtons.emplace_back("AC Sweep Analysis", 100, 100, 180, 30, []() { setup_dialog(DialogType::AC_SWEEP_ANALYSIS, "AC Sweep Analysis", { "Source Name", "Start Freq", "End Freq", "Points/Decade", "Wanted Values" }); isSimulateMenuOpen = false; });
    simulateMenuButtons.emplace_back("Phase Sweep Analysis", 100, 130, 180, 30, []() { setup_dialog(DialogType::PHASE_SWEEP_ANALYSIS, "Phase Sweep Analysis", { "Frequency", "Source Name", "Start Phase", "End Phase", "Increment", "Wanted Values" }); isSimulateMenuOpen = false; });
    vector<Button> resultsViewButtons;
    resultsViewButtons.emplace_back("Back to Editor", 10, SCREEN_HEIGHT - 50, 200, 40, []() { currentAppState = AppState::SCHEMATIC_EDITOR; });
    resultsViewButtons.emplace_back("Probe Voltage", 10, 50, 200, 30, []() { currentInteractionMode = InteractionMode::PROBE_VOLTAGE; });
    resultsViewButtons.emplace_back("Probe Current", 10, 90, 200, 30, []() { currentInteractionMode = InteractionMode::PROBE_CURRENT; });
    resultsViewButtons.emplace_back("Probe Power", 10, 130, 200, 30, []() { currentInteractionMode = InteractionMode::PROBE_POWER; });
    resultsViewButtons.emplace_back("Manage Traces", 10, 170, 200, 30, []() { show_trace_dialog = true; active_trace_dialog_field = 0; SDL_StartTextInput(); });
    resultsViewButtons.emplace_back("Adjust Scale", 10, 210, 200, 30, []() { show_scale_dialog = true; active_trace_dialog_field = 0; SDL_StartTextInput(); });
    vector<ComponentMenuItem> componentMenuItems;
    const int COMP_ITEM_W = 150, COMP_ITEM_H = 100, COMP_ITEM_PAD = 10;
    auto selectComp = [&](ComponentType t) {
        selectedComponentType = t;
        currentInteractionMode = InteractionMode::PLACE_COMPONENT;
        isComponentsMenuOpen = false;
        placementRotation = 0;
    };
    componentMenuItems.push_back({ "Resistor", ComponentType::RESISTOR, loadAndProcessTexture(ASSET_PATH + "resistor.png", renderer), Button("", 210, 40, COMP_ITEM_W, COMP_ITEM_H, [=]() { selectComp(ComponentType::RESISTOR); }) });
    componentMenuItems.push_back({ "Capacitor", ComponentType::CAPACITOR, loadAndProcessTexture(ASSET_PATH + "capacitor.png", renderer), Button("", 210 + (COMP_ITEM_W + COMP_ITEM_PAD), 40, COMP_ITEM_W, COMP_ITEM_H, [=]() { selectComp(ComponentType::CAPACITOR); }) });
    componentMenuItems.push_back({ "Inductor", ComponentType::INDUCTOR, loadAndProcessTexture(ASSET_PATH + "inductor.png", renderer), Button("", 210 + 2 * (COMP_ITEM_W + COMP_ITEM_PAD), 40, COMP_ITEM_W, COMP_ITEM_H, [=]() { selectComp(ComponentType::INDUCTOR); }) });
    componentMenuItems.push_back({ "Diode", ComponentType::DIODE, loadAndProcessTexture(ASSET_PATH + "diode.png", renderer), Button("", 210 + 3 * (COMP_ITEM_W + COMP_ITEM_PAD), 40, COMP_ITEM_W, COMP_ITEM_H, [=]() { selectComp(ComponentType::DIODE); }) });
    componentMenuItems.push_back({ "DC Voltage Src", ComponentType::DC_VOLTAGE_SOURCE, loadAndProcessTexture(ASSET_PATH + "voltage_source.png", renderer), Button("", 210, 40 + (COMP_ITEM_H + COMP_ITEM_PAD), COMP_ITEM_W, COMP_ITEM_H, [=]() { selectComp(ComponentType::DC_VOLTAGE_SOURCE); }) });
    componentMenuItems.push_back({ "DC Current Src", ComponentType::DC_CURRENT_SOURCE, loadAndProcessTexture(ASSET_PATH + "current_source.png", renderer), Button("", 210 + (COMP_ITEM_W + COMP_ITEM_PAD), 40 + (COMP_ITEM_H + COMP_ITEM_PAD), COMP_ITEM_W, COMP_ITEM_H, [=]() { selectComp(ComponentType::DC_CURRENT_SOURCE); }) });
    componentMenuItems.push_back({ "AC Voltage Src", ComponentType::AC_VOLTAGE_SOURCE, loadAndProcessTexture(ASSET_PATH + "ac_voltage_source.png", renderer), Button("", 210 + 2 * (COMP_ITEM_W + COMP_ITEM_PAD), 40 + (COMP_ITEM_H + COMP_ITEM_PAD), COMP_ITEM_W, COMP_ITEM_H, [=]() { selectComp(ComponentType::AC_VOLTAGE_SOURCE); }) });
    auto start_component_placement_from_shortcut = [&](ComponentType t) {
        selectedComponentType = t;
        currentInteractionMode = InteractionMode::PLACE_COMPONENT;
        isComponentsMenuOpen = false;
        isFileMenuOpen = false;
        isSimulateMenuOpen = false;
        placementRotation = 0;
    };
    SDL_Event e;
    while (!quit_flag) {
        while (SDL_PollEvent(&e) != 0) {
            if (e.type == SDL_QUIT) quit_flag = true;
            if (currentAppState == AppState::SCHEMATIC_EDITOR) {
                if (currentInteractionMode == InteractionMode::DIALOG_ACTIVE) {
                    if (e.type == SDL_KEYDOWN) {
                        if (e.key.keysym.sym == SDLK_ESCAPE) close_dialog();
                        if (activeDialogType != DialogType::OPEN) {
                            if (e.key.keysym.sym == SDLK_RETURN) handle_dialog_ok();
                            else if (e.key.keysym.sym == SDLK_TAB && activeDialogFieldIndex != -1) activeDialogFieldIndex = (activeDialogFieldIndex + 1) % dialogFields.size();
                            else if (e.key.keysym.sym == SDLK_BACKSPACE && activeDialogFieldIndex != -1 && !dialogFields[activeDialogFieldIndex].buffer.empty()) dialogFields[activeDialogFieldIndex].buffer.pop_back();
                        }
                    } else if (e.type == SDL_TEXTINPUT && activeDialogType != DialogType::OPEN && activeDialogFieldIndex != -1) {
                        dialogFields[activeDialogFieldIndex].buffer += e.text.text;
                    } else if (e.type == SDL_MOUSEBUTTONDOWN) {
                        int mx, my; SDL_GetMouseState(&mx, &my);
                        if (activeDialogType == DialogType::OPEN) {
                            const int dialog_w = 450, dialog_h = 400;
                            SDL_Rect panelRect = { SCREEN_WIDTH / 2 - dialog_w / 2, SCREEN_HEIGHT / 2 - dialog_h / 2, dialog_w, dialog_h };
                            SDL_Rect cancelRect = { panelRect.x + panelRect.w / 2 - 40, panelRect.y + panelRect.h - 45, 80, 30 };
                            if (SDL_PointInRect(new const SDL_Point{ mx, my }, &cancelRect)) {
                                close_dialog();
                            } else if (dialog_hovered_file_index != -1) {
                                load_schematic(SCHEMATICS_PATH + dialog_file_list[dialog_hovered_file_index]);
                                close_dialog();
                            }
                        } else {
                            for (size_t i = 0; i < dialogFields.size(); ++i) if (SDL_PointInRect(new const SDL_Point{ mx, my }, &dialogFields[i].input_rect)) { activeDialogFieldIndex = i; break; }
                            const int dialog_w = 450, dialog_h = 100 + dialogFields.size() * 40 + 50;
                            SDL_Rect panelRect = { SCREEN_WIDTH / 2 - dialog_w / 2, SCREEN_HEIGHT / 2 - dialog_h / 2, dialog_w, dialog_h };
                            SDL_Rect okRect = { panelRect.x + panelRect.w - 190, panelRect.y + panelRect.h - 45, 80, 30 };
                            SDL_Rect cancelRect = { panelRect.x + panelRect.w - 100, panelRect.y + panelRect.h - 45, 80, 30 };
                            if (SDL_PointInRect(new const SDL_Point{ mx, my }, &okRect)) handle_dialog_ok();
                            else if (SDL_PointInRect(new const SDL_Point{ mx, my }, &cancelRect)) close_dialog();
                        }
                    } else if (e.type == SDL_MOUSEWHEEL && activeDialogType == DialogType::OPEN) {
                        if (e.wheel.y > 0) dialog_file_list_scroll_offset--;
                        else if (e.wheel.y < 0) dialog_file_list_scroll_offset++;
                        if (dialog_file_list_scroll_offset < 0) dialog_file_list_scroll_offset = 0;
                        if (!dialog_file_list.empty() && dialog_file_list_scroll_offset > (int)dialog_file_list.size() - 1) {
                            dialog_file_list_scroll_offset = dialog_file_list.size() - 1;
                        }
                    } else if (e.type == SDL_MOUSEMOTION && activeDialogType == DialogType::OPEN) {
                        int mx, my; SDL_GetMouseState(&mx, &my);
                        const int dialog_w = 450, dialog_h = 400;
                        SDL_Rect panelRect = { SCREEN_WIDTH / 2 - dialog_w / 2, SCREEN_HEIGHT / 2 - dialog_h / 2, dialog_w, dialog_h };
                        int current_y = panelRect.y + 60;
                        const int item_height = 25;
                        const int list_height = dialog_h - 120;
                        const int max_visible_items = list_height / item_height;
                        dialog_hovered_file_index = -1;
                        for (size_t i = 0; i < dialog_file_list.size(); ++i) {
                            int displayed_index = i - dialog_file_list_scroll_offset;
                            if (displayed_index < 0 || displayed_index >= max_visible_items) continue;
                            SDL_Rect item_rect = { panelRect.x + 20, current_y + displayed_index * item_height, panelRect.w - 40, item_height };
                            if (SDL_PointInRect(new const SDL_Point{ mx, my }, &item_rect)) {
                                dialog_hovered_file_index = i;
                                break;
                            }
                        }
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
                    if (e.key.keysym.sym == SDLK_ESCAPE) {
                        currentInteractionMode = InteractionMode::NONE;
                        isDrawingWire = false;
                    }
                    else if (currentInteractionMode == InteractionMode::NONE) {
                        switch(e.key.keysym.sym) {
                            case SDLK_r: start_component_placement_from_shortcut(ComponentType::RESISTOR); break;
                            case SDLK_c: start_component_placement_from_shortcut(ComponentType::CAPACITOR); break;
                            case SDLK_l: start_component_placement_from_shortcut(ComponentType::INDUCTOR); break;
                            case SDLK_d: start_component_placement_from_shortcut(ComponentType::DIODE); break;
                            case SDLK_v: start_component_placement_from_shortcut(ComponentType::DC_VOLTAGE_SOURCE); break;
                            case SDLK_i: start_component_placement_from_shortcut(ComponentType::DC_CURRENT_SOURCE); break;
                            case SDLK_w: currentInteractionMode = InteractionMode::WIRING; break;
                            case SDLK_g: currentInteractionMode = InteractionMode::PLACE_GND_LABEL; break;
                        }
                    }
                    else if (e.key.keysym.sym == SDLK_r && currentInteractionMode == InteractionMode::PLACE_COMPONENT) {
                        placementRotation = (placementRotation + 90) % 360;
                    }
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
                            components.push_back(newComp);
                            schematicModified = true;
                            editingComponent = &components.back();
                            textInputBuffer = editingComponent->value;
                            currentInteractionMode = InteractionMode::EDITING_COMPONENT_VALUE;
                            SDL_StartTextInput();
                        }
                        else if (currentInteractionMode == InteractionMode::DELETE_ITEM) {
                            const double DELETE_THRESHOLD_SQ = 10 * 10;
                            SDL_Point clickPoint = { mouseX, mouseY };
                            int wire_to_delete = -1, component_to_delete = -1, label_to_delete = -1;
                            double min_dist_sq = DELETE_THRESHOLD_SQ;
                            for (int i = 0; i < (int)labels.size(); ++i) if (SDL_PointInRect(&clickPoint, &labels[i].text_rect)) { label_to_delete = i; break; }
                            if (label_to_delete == -1) {
                                for (int i = 0; i < (int)wires.size(); ++i) { double d = dist_to_segment_sq(clickPoint, wires[i].start, wires[i].end); if (d < min_dist_sq) { min_dist_sq = d; wire_to_delete = i; component_to_delete = -1; } }
                                for (int i = 0; i < (int)components.size(); ++i) { double d = dist_to_segment_sq(clickPoint, components[i].node1, components[i].node2); if (d < min_dist_sq) { min_dist_sq = d; component_to_delete = i; wire_to_delete = -1; } }
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
            } else if (currentAppState == AppState::RESULTS_VIEW) {
                if (currentInteractionMode == InteractionMode::PROBE_CURRENT || currentInteractionMode == InteractionMode::PROBE_VOLTAGE || currentInteractionMode == InteractionMode::PROBE_POWER) {
                    if (e.type == SDL_KEYDOWN && e.key.keysym.sym == SDLK_ESCAPE) {
                        currentInteractionMode = InteractionMode::NONE;
                    } else if (e.type == SDL_MOUSEBUTTONDOWN && e.button.button == SDL_BUTTON_LEFT) {
                        int mouseX, mouseY;
                        SDL_GetMouseState(&mouseX, &mouseY);
                        SDL_Point clickPoint = { mouseX, mouseY };
                        string trace_name;
                        if (currentInteractionMode == InteractionMode::PROBE_VOLTAGE) {
                            SDL_Point snapped_click = snap_to_grid(mouseX, mouseY);
                            string ground_node_name;
                            auto point_to_node_name_map = resolve_node_names(ground_node_name);
                            auto point_to_key = [](SDL_Point p) { return to_string(p.x) + "," + to_string(p.y); };
                            string key = point_to_key(snapped_click);
                            if(point_to_node_name_map.count(key)) {
                                trace_name = "V(" + point_to_node_name_map.at(key) + ")";
                            }
                        } else {
                            int component_to_probe = -1;
                            double min_dist_sq = 15 * 15;
                            for (int i = 0; i < (int)components.size(); ++i) {
                                double d = dist_to_segment_sq(clickPoint, components[i].node1, components[i].node2);
                                if (d < min_dist_sq) {
                                    min_dist_sq = d;
                                    component_to_probe = i;
                                }
                            }
                            if (component_to_probe != -1) {
                                if(currentInteractionMode == InteractionMode::PROBE_CURRENT) {
                                    trace_name = "I(" + components[component_to_probe].id + ")";
                                } else {
                                    trace_name = "P(" + components[component_to_probe].id + ")";
                                }
                            }
                        }
                        if (!trace_name.empty()) {
                            bool found = false;
                            for(const auto& pv : plotted_variables) {
                                if (pv.name == trace_name) {
                                    found = true;
                                    break;
                                }
                            }
                            if (!found) {
                                plotted_variables.push_back({trace_name, plot_colors[plotted_variables.size() % plot_colors.size()]});
                            }
                            currentInteractionMode = InteractionMode::NONE;
                        }
                    }
                }
                else if(show_trace_dialog) {
                    if (e.type == SDL_KEYDOWN) {
                        if (e.key.keysym.sym == SDLK_RETURN) {
                            if(!new_trace_buffer.empty()) {
                                plotted_variables.push_back({new_trace_buffer, plot_colors[plotted_variables.size() % plot_colors.size()]});
                                new_trace_buffer = "";
                            }
                        } else if (e.key.keysym.sym == SDLK_ESCAPE) {
                            show_trace_dialog = false;
                            SDL_StopTextInput();
                        } else if (e.key.keysym.sym == SDLK_BACKSPACE && !new_trace_buffer.empty()) {
                            new_trace_buffer.pop_back();
                        }
                    } else if (e.type == SDL_TEXTINPUT) {
                        new_trace_buffer += e.text.text;
                    } else if (e.type == SDL_MOUSEBUTTONDOWN) {
                        int mx, my; SDL_GetMouseState(&mx, &my);
                        const int dialog_w = 400;
                        const int dialog_h = 150 + plotted_variables.size() * 40;
                        SDL_Rect panelRect = { SCREEN_WIDTH / 2 - dialog_w / 2, SCREEN_HEIGHT / 2 - dialog_h / 2, dialog_w, dialog_h };
                        int current_y = panelRect.y + 50;
                        for(size_t i = 0; i < plotted_variables.size(); ++i) {
                            SDL_Rect color_box = {panelRect.x + 20, current_y, 30, 30};
                            if (SDL_PointInRect(new const SDL_Point{mx, my}, &color_box)) {
                                int current_color_index = -1;
                                for(size_t c = 0; c < plot_colors.size(); ++c) {
                                    if(plotted_variables[i].color.r == plot_colors[c].r && plotted_variables[i].color.g == plot_colors[c].g && plotted_variables[i].color.b == plot_colors[c].b) {
                                        current_color_index = c;
                                        break;
                                    }
                                }
                                plotted_variables[i].color = plot_colors[(current_color_index + 1) % plot_colors.size()];
                            }
                            SDL_Rect delete_rect = {panelRect.x + dialog_w - 80, current_y, 60, 30};
                            if (SDL_PointInRect(new const SDL_Point{mx, my}, &delete_rect)) {
                                plotted_variables.erase(plotted_variables.begin() + i);
                                break;
                            }
                            current_y += 40;
                        }
                        SDL_Rect input_rect = {panelRect.x + 60, current_y, dialog_w - 150, 30};
                        if(SDL_PointInRect(new const SDL_Point{mx, my}, &input_rect)) {
                            active_trace_dialog_field = 0;
                        } else {
                            active_trace_dialog_field = -1;
                        }
                        SDL_Rect add_rect = {input_rect.x + input_rect.w + 10, current_y, 60, 30};
                        if (SDL_PointInRect(new const SDL_Point{mx, my}, &add_rect)) {
                            if(!new_trace_buffer.empty()) {
                                plotted_variables.push_back({new_trace_buffer, plot_colors[plotted_variables.size() % plot_colors.size()]});
                                new_trace_buffer = "";
                            }
                        }
                        SDL_Rect close_rect = {panelRect.x + dialog_w/2 - 40, panelRect.y + dialog_h - 45, 80, 30};
                        if (SDL_PointInRect(new const SDL_Point{mx, my}, &close_rect)) {
                            show_trace_dialog = false;
                            SDL_StopTextInput();
                        }
                    }
                } else if (show_scale_dialog) {
                    vector<tuple<string, string*, int>> fields = {
                            make_tuple("X-Min:", &manual_scale.x_min_str, 0),
                            make_tuple("X-Max:", &manual_scale.x_max_str, 1),
                            make_tuple("Y-Min:", &manual_scale.y_min_str, 2),
                            make_tuple("Y-Max:", &manual_scale.y_max_str, 3)
                    };
                    if (e.type == SDL_KEYDOWN) {
                        if (e.key.keysym.sym == SDLK_RETURN || e.key.keysym.sym == SDLK_ESCAPE) {
                            show_scale_dialog = false;
                            SDL_StopTextInput();
                        } else if (e.key.keysym.sym == SDLK_TAB) {
                            active_trace_dialog_field = (active_trace_dialog_field + 1) % fields.size();
                        } else if (e.key.keysym.sym == SDLK_BACKSPACE && active_trace_dialog_field != -1) {
                            string* buffer = get<1>(fields[active_trace_dialog_field]);
                            if (!buffer->empty()) buffer->pop_back();
                        }
                    } else if (e.type == SDL_TEXTINPUT && active_trace_dialog_field != -1) {
                        *get<1>(fields[active_trace_dialog_field]) += e.text.text;
                    } else if (e.type == SDL_MOUSEBUTTONDOWN) {
                        int mx, my; SDL_GetMouseState(&mx, &my);
                        const int dialog_w = 350;
                        const int dialog_h = 280;
                        SDL_Rect panelRect = { SCREEN_WIDTH / 2 - dialog_w / 2, SCREEN_HEIGHT / 2 - dialog_h / 2, dialog_w, dialog_h };
                        int current_y = panelRect.y + 50;
                        active_trace_dialog_field = -1;
                        for(const auto& field : fields) {
                            SDL_Rect input_rect = {panelRect.x + 80, current_y, dialog_w - 100, 30};
                            if(SDL_PointInRect(new const SDL_Point{mx, my}, &input_rect)) {
                                active_trace_dialog_field = get<2>(field);
                            }
                            current_y += 40;
                        }
                        SDL_Rect close_rect = {panelRect.x + dialog_w/2 - 40, panelRect.y + dialog_h - 45, 80, 30};
                        if (SDL_PointInRect(new const SDL_Point{mx, my}, &close_rect)) {
                            show_scale_dialog = false;
                            SDL_StopTextInput();
                        }
                    }
                } else {
                    if (e.type == SDL_KEYDOWN) {
                        if (e.key.keysym.sym == SDLK_RIGHT) {
                            if (last_simulated_circuit) {
                                const bool is_tran = last_simulated_circuit->tran_solved;
                                const bool is_dc = last_simulated_circuit->dc_sweep_solved;
                                const bool is_ac = last_simulated_circuit->ac_sweep_solved;
                                const bool is_phase = last_simulated_circuit->phase_sweep_solved;
                                const auto& results = is_tran ? last_simulated_circuit->tran_results :
                                                      is_dc ? last_simulated_circuit->dc_sweep_results :
                                                      is_ac ? last_simulated_circuit->ac_sweep_results :
                                                      last_simulated_circuit->phase_sweep_results;
                                if (!results.empty() && cursor_index < (int)results.size() - 1) {
                                    cursor_index++;
                                }
                            }
                        } else if (e.key.keysym.sym == SDLK_LEFT) {
                            if (cursor_index > 0) {
                                cursor_index--;
                            }
                        }
                    }
                    for (auto& b : resultsViewButtons) {
                        b.handle_event(&e);
                    }
                }
            }
        }
        SDL_SetRenderDrawColor(renderer, 0x22, 0x22, 0x22, 0xFF);
        SDL_RenderClear(renderer);
        if (currentAppState == AppState::SCHEMATIC_EDITOR) {
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
        } else if (currentAppState == AppState::RESULTS_VIEW) {
            render_results_view(renderer, uiFont, resultsViewButtons, componentMenuItems);
        }
        SDL_RenderPresent(renderer);
    }
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