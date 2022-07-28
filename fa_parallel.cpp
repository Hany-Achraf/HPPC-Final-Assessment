/*
* Hany Ashraf Mohamed (A18CS4012) - Sec 02
* HPPC Final Assessment Source Code (Part 1)
*/

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <mpi.h>
using namespace std;

// Square struct that basically will be used to hold the four coordinates of a square to ease handling squares
struct Square {
  // Attributes (the four coordinates)
  int top_left_coordinate1;
  int top_left_coordinate2;
  int bottom_right_coordinate1;
  int bottom_right_coordinate2;
  // Constructor
  Square(int top_left_coordinate1, int top_left_coordinate2, int bottom_right_coordinate1, int bottom_right_coordinate2) {
        this->top_left_coordinate1 = top_left_coordinate1;
        this->top_left_coordinate2 = top_left_coordinate2;
        this->bottom_right_coordinate1 = bottom_right_coordinate1;
        this->bottom_right_coordinate2 = bottom_right_coordinate2;
  }
};

// This func takes two squares as inputs and returns true if the first square has the second one as a sub-matrix
bool has_sub_square(Square *sqr1, Square *sqr2) {
      if (sqr2->top_left_coordinate1 >= sqr1->top_left_coordinate1 &&
          sqr2->top_left_coordinate2 >= sqr1->top_left_coordinate2 &&
          sqr2->bottom_right_coordinate1 <= sqr1->bottom_right_coordinate1 &&
          sqr2->bottom_right_coordinate2 <= sqr1->bottom_right_coordinate2) return true;
      else
          return false;
}

// This func takes a 1D array and convert it to 2D array (matrix)
void convert_1D_arr_to_2D_matrix(int num_elemts_per_proc, char arr[], char matrix[][5000]) {
  for (int i = 0; i < num_elemts_per_proc; i++)
    matrix[i / 5000][i % 5000] = arr[i];
}

// This func finds the squares within a sub-matrix, and the row index to be sent another process for further processing.
// (Output parameters: sqrs, row_index)
void process_matrix(char sub_matrix[][5000], int elements_per_proc, int num_recv_elemts, int proc_rank,
                    int last_working_proc_rank, vector<Square> &sqrs, int* row_index) {
  bool out_of_bound = false, continue_loop;
  int x, y, x_count, y_count;

  for (int i = 0; i < (elements_per_proc + num_recv_elemts) / 5000 && !out_of_bound; i++) {
    for (int j = 0; j < 5000 && !out_of_bound; j++) {
      if (sub_matrix[i][j] == '2') {
          x = i + 1, y = j + 1, continue_loop = false;
          if (x < 5000 && y < 5000)
            continue_loop = true;
          while (x < (elements_per_proc + num_recv_elemts) / 5000 && y < 5000 && continue_loop) {
            for (x_count = i; x_count <= x && continue_loop; x_count++)
                for (y_count = j; y_count <= y && continue_loop; y_count++)
                  if (sub_matrix[x_count][y_count] != '2')
                    continue_loop = false;

            if (continue_loop && x == (elements_per_proc + num_recv_elemts) / 5000 - 1 && proc_rank == last_working_proc_rank) {
              continue_loop = false;
              x++, y++;
            }

            if (!continue_loop) {
              x--, y--;
              if (x - i == y - j && x - i > 0) {
                bool not_sub_sqr = true;

                Square *sqr;
                if (num_recv_elemts == 0)
                  sqr = new Square(i + proc_rank * (elements_per_proc / 5000), j,
                                   x + proc_rank * (elements_per_proc / 5000), y);
                else
                  sqr = new Square(i + (elements_per_proc / 5000 * proc_rank) - (num_recv_elemts / 5000), j,
                                   x + (elements_per_proc / 5000 * proc_rank) - (num_recv_elemts / 5000), y);

                for (int count = 0; count < sqrs.size(); count++)
                    if (has_sub_square(&sqrs[count], sqr))
                        not_sub_sqr = false;
                if (not_sub_sqr)
                    sqrs.push_back(*sqr);
              }
            }
            x++, y++;
          }
          if (continue_loop) {
            *row_index = i;
            out_of_bound = true;
          }
      }
    }
  }
}

// This func takes a filename as an input, and fill arr (output argument) with the elements in the file
void read_data_from_file(string filename, char arr[]) {
  ifstream input_file(filename);
  if (!input_file.is_open()) {
      cerr << "Could not open the file - '" << filename << "'" << endl;
      return;
  }

  int i = 0;
  while (!input_file.eof())
      input_file >> arr[i++];

  input_file.close();
}

// this array will hold all the read data from the file (Defined outside of the main func, because its length is more than 1000000)
char arr[25000000];
// this array will be used later to receive elements during P2P communication among the processes
char recv_buf_from_proc[5000000];

int main(int argc, char** argv) {
  // Filling arr with the elements in the file "data.txt"
  // MUST BE IN THE SAME FOLDER DURING RUNNING, OR YOU MAY CHANGE THE PATH STRING)
  read_data_from_file("data.txt", arr);

  int myrank, nps, elements_per_proc, num_recv_elemts = 0;

  // Starting MPI parallelization
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &nps);
  MPI_Status status;
  MPI_Request request;

  //Defining a new MPI datatype for Square struct
  MPI_Datatype dt_square;
  MPI_Type_contiguous(4, MPI_INT, &dt_square);
  MPI_Type_commit(&dt_square);

  elements_per_proc = 25000000 / nps;
  while (elements_per_proc % 5000 != 0)
    elements_per_proc++;

  // last_working_proc_rank equals nps - 1, if nps is a divisor for 5000. Otherwise it equals the rank of the last process that
  // will receive elements from the 25000000, since in this case there will some be processes at the end will hold some
  // data outside that don't belong to the array (Note that we increase the number of elements_per_proc till its a divisor of 5000)
  int last_working_proc_rank = ceil((float)25000000 / elements_per_proc) - 1;

  // Create a buffer that will hold a subset of the numbers
  char* recv_buf = (char*)malloc(sizeof(int) * elements_per_proc);

  // Scatter the numbers to all processes
  MPI_Scatter(arr, elements_per_proc, MPI_CHAR, recv_buf, elements_per_proc, MPI_CHAR, 0, MPI_COMM_WORLD);

  // Checking whetherthe process has some data from the 25000000 to be processed or not. (If nps 10, last_working_proc_rank = nps - 1)
  if (myrank <= last_working_proc_rank) {
    // Defining a buffer called new_buf, that will be used later in case my process wants to communicate with a neighbor process
    // Initializing it with the data in the recv_buf
    char* new_buf = (char*)malloc(sizeof(char) * elements_per_proc);
    for (int i = 0; i < elements_per_proc; i++)
      new_buf[i] = recv_buf[i];

    // This vector will hold my squares, as a process, that I will send it back to the master at the end of pararllelization
    vector<Square> sqrs;
    // This vector will hold the edge squares discovered that I will send it back in response to a neighbor process
    vector<Square> edge_sqrs;

    // In case that I reach to an element that needs more rows (processing by a neighbor process),
    // I need to know the row index of this element, since all all the elements starting from it will be sent to the neighbor process for processing
    int row_index = -1;
    int recv_row_index = -1;

    // All the processes, except the master, will receive a row index from the process with rank less than by 1
    if (myrank > 0)
      MPI_Irecv(&recv_row_index, 1, MPI_INT, myrank - 1, 0, MPI_COMM_WORLD, &request);

    // Creating and filling a sub-matrix using the elements in the recevied buffer
    char sub_matrix[elements_per_proc / 5000][5000];
    convert_1D_arr_to_2D_matrix(elements_per_proc, recv_buf, sub_matrix);

    // Processing this sub matrix (You may refer the comment above this func)
    process_matrix(sub_matrix, elements_per_proc, 0, myrank, last_working_proc_rank, sqrs, &row_index);

    // Waiting the row index coming from the process with rank = myrank - 1.
    if (myrank > 0) {
      MPI_Wait(&request, &status);

      // and then start receiving the actual elements that exist in this row and further to be processed by me, as a process.
      if (recv_row_index > -1) {
        MPI_Recv(recv_buf_from_proc, elements_per_proc, MPI_CHAR, myrank - 1, 1, MPI_COMM_WORLD, &status);
        MPI_Get_count(&status, MPI_CHAR, &num_recv_elemts);

        // Updating my row index after receiving the elements, because it may happen that I also need to pass these elements or part of it
        // plus my elements to the next process
        if (row_index > -1)
          row_index += (num_recv_elemts / 5000);

        // Initializing and filling a new buffer rather than the recevied one, the new buffer will contain the elements I have just received
        // plus my elements that I initially recvied on MPI_Scatter()
        new_buf = (char*)malloc(sizeof(char) * (elements_per_proc + num_recv_elemts));
        for (int i = 0; i < num_recv_elemts; i++)
          new_buf[i] = recv_buf_from_proc[i];
        for (int i = 0; i < elements_per_proc; i++)
            new_buf[i + num_recv_elemts] = recv_buf[i];

        // Creating and filling a new sub-matrix using the elements in the new buffer
        char new_matrix[(elements_per_proc + num_recv_elemts) / 5000][5000];
        convert_1D_arr_to_2D_matrix(elements_per_proc + num_recv_elemts, new_buf, new_matrix);

        // Processing this new matrix (You may refer the comment associated with this func)
        process_matrix(new_matrix, elements_per_proc, num_recv_elemts, myrank, last_working_proc_rank, edge_sqrs, &row_index);

      }
    }

    // If my rank is less than last working process rank (which normally equals nps unless nps is NOT a divisor for 5000),
    // send the row index that has elements need more processing in process with a rank = myrank + 1
    if (myrank < last_working_proc_rank) {
      MPI_Send(&row_index, 1, MPI_INT, myrank + 1, 0, MPI_COMM_WORLD);
      if (row_index > -1) {
        // Sending the elements starting from this row index in my sub-marix to the process (myrank + 1)
        // Note that we send the new buffer because a square might be shredded among more than 2 sub-marices/processes
        MPI_Send(new_buf + (row_index * 5000), elements_per_proc + num_recv_elemts - row_index * 5000, MPI_CHAR, myrank + 1, 1, MPI_COMM_WORLD);

        // Receving the edge squares coming in response to me sending some of my elements to another process
        Square *recv_sqrs_from_proc = (Square*)malloc(sizeof(Square) * 100); // Assuming that the maxmum num of recevied squares is less than 100
        MPI_Recv(recv_sqrs_from_proc, 100, dt_square, myrank + 1, 2, MPI_COMM_WORLD, &status);

        // Getting the acual number of the squares just received
        int num_recv_sqrs;
        MPI_Get_count(&status, dt_square, &num_recv_sqrs);

        for (int i = 0; i < num_recv_sqrs; i++)
          // Checking if the recevied edge square belongs to me or not, if NO, it will be pushed to the edge_sqrs vector
          if (recv_sqrs_from_proc[i].top_left_coordinate1 < (elements_per_proc / 5000) * myrank) {
            for (int j = 0; j < edge_sqrs.size(); j++)
              if (has_sub_square(&recv_sqrs_from_proc[i], &edge_sqrs[j]))
                edge_sqrs.erase(edge_sqrs.begin() + j--);
            edge_sqrs.push_back(recv_sqrs_from_proc[i]);
          } else {
            // If the recevied edge square belongs to me, so it will be pushed to my sqrs.
            sqrs.push_back(recv_sqrs_from_proc[i]);
          }
      }
    }


    if (recv_row_index > -1) {
      // In case I received elements from a process with a rank less by 1
      // If I found edge squares that don't start/belong to me, I will delete all my squares that are sub-squares within the edge squares
      for (int i = 0; i < edge_sqrs.size(); i++)
        for (int j = 0; j < sqrs.size(); j++)
          if (has_sub_square(&edge_sqrs[i], &sqrs[j]))
            sqrs.erase(sqrs.begin() + j--);

      // Then I will send the edge squares to the process that they belong to
      MPI_Send(&edge_sqrs[0], edge_sqrs.size(), dt_square, myrank - 1, 2, MPI_COMM_WORLD);
    }

    // Finally the master process will receive all the squares from the processes, and present them.
    if (myrank == 0) {
      Square* all_sqrs_buf = (Square*)malloc(sizeof(Square) * 1000);
      vector<Square> all_sqrs;

      for (int i = 0; i < sqrs.size(); i++)
        all_sqrs.push_back(sqrs[i]);

      int num_recv_sqrs = 0;
      for (int i = 1; i <= last_working_proc_rank; i++) {
        MPI_Recv(all_sqrs_buf, 1000, dt_square, i, 5, MPI_COMM_WORLD, &status);
        MPI_Get_count(&status, dt_square, &num_recv_sqrs);
        for (int j = 0; j < num_recv_sqrs; j++)
          all_sqrs.push_back(all_sqrs_buf[j]);
      }

      cout << "Number of squares: " << all_sqrs.size() << "\n\n";
      for (int i = 0; i < all_sqrs.size(); i++)
      cout << "Square " << i + 1 << " coordinates [top_left_corner]: " << "("
                                 << all_sqrs[i].top_left_coordinate1 << "," << all_sqrs[i].top_left_coordinate2 << ")" << "\n"
           << "Square " << i + 1 << " coordinates [bottom_right_corner]: " << "("
                                 << all_sqrs[i].bottom_right_coordinate1 << "," << all_sqrs[i].bottom_right_coordinate2 << ")\n\n";
    } else {
      // Sending all the discovered squares, as a one process, to the master process for counting and presenting all
      MPI_Send(&sqrs[0], sqrs.size(), dt_square, 0, 5, MPI_COMM_WORLD);
    }
  }

  MPI_Finalize();
  return 0;
}
