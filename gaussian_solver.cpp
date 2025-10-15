#include <bits/stdc++.h>
using namespace std;

static constexpr double EPS = 1e-10;

struct Logger {
    int iter = 1;
    int n = 0; // rows (square system)
    int cols = 0; // n + 1 (augmented)
    bool show_row_before_after = true;

    void print_row(const vector<double>& row, int cols) {
        cout << "[ ";
        cout << fixed << setprecision(6);
        for (int j = 0; j < cols - 1; ++j) {
            double v = (fabs(row[j]) < 1e-12 ? 0.0 : row[j]);
            cout << setw(12) << v << " ";
        }
        double b = (fabs(row[cols-1]) < 1e-12 ? 0.0 : row[cols-1]);
        cout << "| " << setw(12) << b << " ]";
    }

    void print_aug(const vector<vector<double>>& M) {
        cout << "Matrix [A | b]:\n";
        for (int i = 0; i < n; ++i) {
            cout << "  ";
            print_row(M[i], cols);
            cout << "\n";
        }
        cout << "\n";
    }

    void op_header(const string& what) {
        cout << "Iter " << (iter++) << ": " << what << "\n";
    }

    void show_before_after_row(int r, const vector<double>& before, const vector<double>& after) {
        if (!show_row_before_after) return;
        cout << "    Row " << (r+1) << " before: ";
        print_row(before, cols);
        cout << "\n";
        cout << "    Row " << (r+1) << " after : ";
        print_row(after, cols);
        cout << "\n";
    }
};

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int n;
    cout << "Enter n (size of A is n x n): ";
    if (!(cin >> n) || n <= 0) {
        cerr << "Invalid n.\n";
        return 1;
    }

    vector<vector<double>> A(n, vector<double>(n + 1));
    cout << "Enter A (n x n) row-wise, then b (n x 1). Total "
         << n << " rows, each with " << (n + 1) << " values:\n";
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j <= n; ++j) {
            if (!(cin >> A[i][j])) {
                cerr << "Invalid input.\n";
                return 1;
            }
        }
    }

    Logger log;
    log.n = n;
    log.cols = n + 1;

    cout << "\nInitial matrix:\n";
    log.print_aug(A);

    for (int col = 0; col < n; ++col) {
        // 1) Partial pivot: find row with max |A[r][col]| for r >= col
        int pivot = col;
        double best = fabs(A[pivot][col]);
        for (int r = col + 1; r < n; ++r) {
            double v = fabs(A[r][col]);
            if (v > best) best = v, pivot = r;
        }

        if (best < EPS) {
            log.op_header("No usable pivot in column " + to_string(col+1) +
                          " (≈ 0). System may be singular or underdetermined. Stopping.");
            log.print_aug(A);
            return 0;
        }

        // 2) Swap pivot row into place (if needed)
        if (pivot != col) {
            log.op_header("Swap R" + to_string(col+1) + " ↔ R" + to_string(pivot+1));
            swap(A[pivot], A[col]);
            log.print_aug(A);
        }

        // 3) Scale pivot row to make pivot = 1
        {
            double p = A[col][col];
            if (fabs(p) < EPS) {
                log.op_header("Unexpected near-zero pivot after swap. Aborting.");
                return 0;
            }
            double inv = 1.0 / p;

            vector<double> before = A[col];
            for (int j = col; j <= n; ++j) A[col][j] *= inv;

            log.op_header("Scale R" + to_string(col+1) + " by " + to_string(inv) +
                          " to make pivot 1 (pivot at row " + to_string(col+1) +
                          ", col " + to_string(col+1) + ")");
            log.show_before_after_row(col, before, A[col]);
            log.print_aug(A);
        }

        // 4) Eliminate column in all other rows
        for (int r = 0; r < n; ++r) {
            if (r == col) continue;
            double factor = A[r][col];
            if (fabs(factor) < EPS) continue;

            vector<double> before = A[r];
            for (int j = col; j <= n; ++j) {
                A[r][j] -= factor * A[col][j];
            }

            // Note: canonical elementary op is Rr ← Rr - factor * Rcol
            log.op_header("R" + to_string(r+1) + " ← R" + to_string(r+1) +
                          " - (" + to_string(factor) + ") * R" + to_string(col+1));
            log.show_before_after_row(r, before, A[r]);
            log.print_aug(A);
        }
    }

    // After Gauss–Jordan, left block ~ I, right is the solution x
    cout << "Finished Gauss–Jordan. Final matrix should be [I | x]:\n";
    log.print_aug(A);

    cout << "Solution vector x:\n";
    cout << fixed << setprecision(10);
    for (int i = 0; i < n; ++i) {
        double v = (fabs(A[i][n]) < 1e-12 ? 0.0 : A[i][n]);
        cout << "x[" << i << "] = " << v << "\n";
    }

    return 0;
}
