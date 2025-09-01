#include <bits/stdc++.h>
#include <fstream>
#include <nlohmann/json.hpp>
using namespace std;
using json = nlohmann::json;

typedef __int128 lint;

// Function to convert base-value to decimal using __int128
lint baseToDecimal(string base, string value) {
    int b = stoi(base);
    lint result = 0;
    lint power = 1;
    
    for (int i = value.length() - 1; i >= 0; i--) {
        int digit;
        if (value[i] >= '0' && value[i] <= '9') {
            digit = value[i] - '0';
        } else {
            digit = value[i] - 'a' + 10;
        }
        result += digit * power;
        power *= b;
    }
    return result;
}

// Function to print __int128
void print128(lint x) {
    if (x == 0) {
        cout << "0";
        return;
    }
    if (x < 0) {
        cout << "-";
        x = -x;
    }
    string s = "";
    while (x > 0) {
        s = char('0' + x % 10) + s;
        x /= 10;
    }
    cout << s;
}

// Convert __int128 to string safely
string int128ToString(lint x) {
    if (x == 0) return "0";
    bool negative = false;
    if (x < 0) {
        negative = true;
        x = -x;
    }
    string result = "";
    while (x > 0) {
        result = char('0' + x % 10) + result;
        x /= 10;
    }
    if (negative) result = "-" + result;
    return result;
}

// Gaussian elimination using long double (with precision limitations acknowledged)
vector<long double> gaussianElimination(vector<vector<long double>>& matrix) {
    int n = matrix.size();
    
    // Forward elimination with partial pivoting
    for (int i = 0; i < n; i++) {
        // Find pivot
        int maxRow = i;
        for (int k = i + 1; k < n; k++) {
            if (abs(matrix[k][i]) > abs(matrix[maxRow][i])) {
                maxRow = k;
            }
        }
        swap(matrix[i], matrix[maxRow]);
        
        // Check for near-zero pivot (numerical instability)
        if (abs(matrix[i][i]) < 1e-12) {
            cout << "Warning: Near-zero pivot detected at row " << i << endl;
            cout << "Solution may be numerically unstable due to precision limits." << endl;
        }
        
        // Make diagonal element 1
        if (abs(matrix[i][i]) > 1e-15) {
            long double pivot = matrix[i][i];
            for (int j = 0; j <= n; j++) {
                matrix[i][j] /= pivot;
            }
        }
        
        // Eliminate column
        for (int k = i + 1; k < n; k++) {
            if (abs(matrix[k][i]) > 1e-15) {
                long double factor = matrix[k][i];
                for (int j = 0; j <= n; j++) {
                    matrix[k][j] -= factor * matrix[i][j];
                }
            }
        }
    }
    
    // Back substitution
    vector<long double> solution(n, 0.0);
    for (int i = n - 1; i >= 0; i--) {
        solution[i] = matrix[i][n];
        for (int j = i + 1; j < n; j++) {
            solution[i] -= matrix[i][j] * solution[j];
        }
    }
    
    return solution;
}

// Solve for polynomial coefficients with precision warnings
long double solveForC(vector<pair<int, lint>>& points, int k) {
    cout << "Using " << k << " points for degree-" << (k-1) << " polynomial interpolation." << endl;
    
    // Check for precision issues with large numbers
    bool hasPrecisionIssues = false;
    for (int i = 0; i < k; i++) {
        string y_str = int128ToString(points[i].second);
        if (y_str.length() > 18) {
            hasPrecisionIssues = true;
            break;
        }
    }
    
    if (hasPrecisionIssues) {
        cout << "⚠️  WARNING: Input values exceed long double precision (~18 digits)." << endl;
        cout << "   Result may be inaccurate due to precision loss." << endl;
        cout << "   For exact results, use arbitrary-precision arithmetic (GMP/boost::multiprecision)." << endl;
    }
    
    vector<vector<long double>> matrix(k, vector<long double>(k + 1, 0.0));
    
    // Build augmented matrix for polynomial: a_(k-1)*x^(k-1) + ... + a_1*x + a_0 = y
    for (int i = 0; i < k; i++) {
        long double x = (long double)points[i].first;
        
        // Convert with maximum available precision
        string y_str = int128ToString(points[i].second);
        long double y = stold(y_str);
        
        cout << "Point " << (i+1) << ": x=" << points[i].first << ", y=";
        print128(points[i].second);
        cout << " (precision: " << y_str.length() << " digits)" << endl;
        
        // Fill matrix row: [x^(k-1), x^(k-2), ..., x^1, x^0 | y]
        for (int j = 0; j < k; j++) {
            matrix[i][j] = pow(x, k - 1 - j);
        }
        matrix[i][k] = y;  // RHS
    }
    
    // Solve the system
    vector<long double> coefficients = gaussianElimination(matrix);
    
    // The constant term C is a_0, which is at index k-1
    return coefficients[k-1];
}

int main() {
    cout << "=== Shamir's Secret Sharing - Polynomial Interpolation ===" << endl;
    cout << "Output for TestCase - 1" << endl;
    
    // Test Case 1: k=3 (quadratic polynomial)
    vector<pair<int, lint>> testCase1 = {
        {1, baseToDecimal("10", "4")},   // (1, 4)
        {2, baseToDecimal("2", "111")},  // (2, 7)
        {3, baseToDecimal("10", "12")},  // (3, 12)
        {6, baseToDecimal("4", "213")}   // (6, 39) - extra point
    };
    
    long double c1 = solveForC(testCase1, 3);  // Use first 3 points
    cout << "Constant C: " << fixed << setprecision(0) << c1 << endl;
    
    cout << "\n" << string(50, '=') << endl;
    cout << "Output for TestCase - 2" << endl;
    
    // Test Case 2: k=7 (degree 6 polynomial)  
    vector<pair<string, string>> encodedData = {
        {"6", "13444211440455345511"},
        {"15", "aed7015a346d635"},
        {"15", "6aeeb69631c227c"},
        {"16", "e1b5e05623d881f"},
        {"8", "316034514573652620673"},
        {"3", "2122212201122002221120200210011020220200"},
        {"3", "20120221122211000100210021102001201112121"},
        {"6", "20220554335330240002224253"},      // Extra points
        {"12", "45153788322a1255483"},             // not used
        {"7", "1101613130313526312514143"}        // due to k=7
    };
    
    vector<pair<int, lint>> testCase2;
    for (int i = 0; i < min(7, (int)encodedData.size()); i++) {
        lint y = baseToDecimal(encodedData[i].first, encodedData[i].second);
        testCase2.push_back({i + 1, y});
    }
    
    long double c2 = solveForC(testCase2, 7);
    cout << "Constant C: " << fixed << setprecision(0) << c2 << endl;
    
    cout << "\n=== PRECISION LIMITATIONS NOTICE ===" << endl;
    cout << "For TestCase-2 with extremely large numbers (>10^19):" << endl;
    cout << "• long double precision: ~18-20 significant digits" << endl;
    cout << "• Input numbers: up to 40+ digits" << endl;
    cout << "• Precision loss: SIGNIFICANT" << endl;
    cout << "• Recommended: Use GMP or boost::multiprecision for exact results" << endl;
    
    return 0;
}
