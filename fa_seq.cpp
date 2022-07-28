#include <iostream>
#include <fstream>
#include <vector>
using namespace std;

struct Square {
  int top_left_coordinate1;
  int top_left_coordinate2;
  int bottom_right_coordinate1;
  int bottom_right_coordinate2;
  Square(int top_left_coordinate1, int top_left_coordinate2, int bottom_right_coordinate1, int bottom_right_coordinate2) {
        this->top_left_coordinate1 = top_left_coordinate1;
        this->top_left_coordinate2 = top_left_coordinate2;
        this->bottom_right_coordinate1 = bottom_right_coordinate1;
        this->bottom_right_coordinate2 = bottom_right_coordinate2;
  }
};

bool has_sub_square(Square *sqr1, Square *sqr2) {
      if (sqr2->top_left_coordinate1 >= sqr1->top_left_coordinate1 &&
          sqr2->top_left_coordinate2 >= sqr1->top_left_coordinate2 &&
          sqr2->bottom_right_coordinate1 <= sqr1->bottom_right_coordinate1 &&
          sqr2->bottom_right_coordinate2 <= sqr1->bottom_right_coordinate2) return true;
      else
          return false;
}

char a[25000000];
char arr[5000][5000];

int main() {
    ifstream input_file("data.txt");

    char ch = 0;
    int i = 0;
    while (!input_file.eof()) {
        input_file >> ch;
        a[i] = ch;
        i++;
    }

    input_file.close();

    for (int i = 0; i < 25000000; i++)
      arr[i / 5000][i % 5000] = a[i];
      
    vector<Square> sqrs;


    bool continue_loop;
    int x, y, sqrs_count = 0;

    for (int i = 0; i < 5000; i++) { // decrease by 1
        for (int j = 0; j < 5000; j++) {  // decrease by 1
            if (arr[i][j] == '2') {
                x = i, y = j, continue_loop = true;
                while (x < 5000 && y < 5000 && continue_loop) {
                    x++, y++;
                    for (int x_count = i; x_count <= x && continue_loop; x_count++)
                        for (int y_count = j; y_count <= y && continue_loop; y_count++)
                            if (arr[x_count][y_count] != '2')
                                continue_loop = false;

                }
                x--, y--;
                if (x - i == y - j && x - i > 0) {
                    bool not_sub_sqr = true;

                    Square sqr(i, j, x, y);
                    for (int count = 0; count < sqrs.size(); count++)
                        if (has_sub_square(&sqrs[count], &sqr))
                            not_sub_sqr = false;

                    if (not_sub_sqr) {
                        Square sqr(i, j, x, y);
                        sqrs.push_back(sqr);
                        cout << "Square #" << ++sqrs_count << " => "
                             << "[top_left_corner]: " << "(" << sqr.top_left_coordinate1 << ","  << sqr.top_left_coordinate2 << ")" << " - "
                             << "[bottom_right_corner]: " << "(" << sqr.bottom_right_coordinate1 << "," << sqr.bottom_right_coordinate2 << ")" << endl;
                    }
                }
            }
        }
    }

    return 0;
}
