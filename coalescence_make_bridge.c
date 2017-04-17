#include<stdio.h>
#include<stdlib.h>
#include<assert.h>
#include<math.h>

#define NUMBER_OF_CRACKS 10
#define DIMENSION 3
#define NUMBER_OF_CRACKFRONT_POINTS 10000

int coalescence_count;
int number_of_coalescence_points;
int coalescence_crack[NUMBER_OF_CRACKS][2];
int cracks;
int nnodes[NUMBER_OF_CRACKS];
double nodes_coordinate[NUMBER_OF_CRACKS][NUMBER_OF_CRACKFRONT_POINTS][DIMENSION];
double coalesced_nodes_coordinate[NUMBER_OF_CRACKFRONT_POINTS][DIMENSION];
double nodal_distance[NUMBER_OF_CRACKS];
double init_univec_normal[NUMBER_OF_CRACKS][DIMENSION];
double p_to_p_vector[NUMBER_OF_CRACKS][NUMBER_OF_CRACKFRONT_POINTS][DIMENSION];
double univec_normal[NUMBER_OF_CRACKS][NUMBER_OF_CRACKFRONT_POINTS][DIMENSION];
double univec_tangent[NUMBER_OF_CRACKS][NUMBER_OF_CRACKFRONT_POINTS][DIMENSION];
double univec_propa[NUMBER_OF_CRACKS][NUMBER_OF_CRACKFRONT_POINTS][DIMENSION];
int temp_number_of_coalescence_points[NUMBER_OF_CRACKS];
double temp_coalescence_points_coordinte[NUMBER_OF_CRACKS][NUMBER_OF_CRACKFRONT_POINTS][DIMENSION];




double Distance(double a[], double b[])
{
  double dist;
  //printf("a = %lf %lf %lf\nb = %lf %lf %lf\n", a[0], a[1], a[2], b[0], b[1], b[2]);
  dist = sqrt((a[0]-b[0])*(a[0]-b[0])+(a[1]-b[1])*(a[1]-b[1])+(a[2]-b[2])*(a[2]-b[2]));
  return dist;
}

double InnerProduct(double a[], double b[])
{
  double product;
  product = a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
  return product;
}

/*
 * ベクトルa,bの外積ベクトルcを算出
 * (a×b=c)
 */
void CrossProduct(double a[], double b[], double c[])
{
  c[0] = a[1]*b[2] - a[2]*b[1];
  c[1] = a[2]*b[0] - a[0]*b[2];
  c[2] = a[0]*b[1] - a[1]*b[0];
}

/*
 * 入力したベクトルを単位ベクトルへと変換する
 */
void GetUnitVector(double a[])
{
  double temp_amount;
  //printf("##before a = (%lf, %lf, %lf)\n", a[0], a[1], a[2]);
  temp_amount = sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
  //printf("##temp_amount = %lf\n", temp_amount);

  a[0] = a[0]/temp_amount;
  a[1] = a[1]/temp_amount;
  a[2] = a[2]/temp_amount;
  //printf("##after a = (%lf, %lf, %lf)\n", a[0], a[1], a[2]);
}

void ReadCoalescenceFlag(const char *filename)
{
  FILE *fp;

  int i;

  fp = fopen(filename, "r");
  assert(fp != NULL);

  fscanf(fp, "%d", &coalescence_count);

  for(i = 0; i < coalescence_count; i++){
    fscanf(fp, "%d %d", &coalescence_crack[i][0], &coalescence_crack[i][1]);
  }
  fclose(fp);
}

void ReadCrackFrontPoints(const char *filename, int crack_count)
{
  FILE *fp;

  int i;
  int dummy;

  fp = fopen(filename, "r");
  assert(fp != NULL);

  fscanf(fp, "%d", &nnodes[crack_count]);

  for(i = 0; i < nnodes[crack_count]; i++){
    fscanf(fp, "%d %lf %lf %lf", &dummy, &nodes_coordinate[crack_count][i][0], &nodes_coordinate[crack_count][i][1], &nodes_coordinate[crack_count][i][2]);
  }
  fclose(fp);
}

void ReadCrackParam(const char *filename, int crack_count)
{
  FILE *fp;

  fp = fopen(filename, "r");
  assert(fp != NULL);

  fscanf(fp, "%lf", &nodal_distance[crack_count]);
  fscanf(fp, "%lf %lf %lf", &init_univec_normal[crack_count][0], &init_univec_normal[crack_count][1], &init_univec_normal[crack_count][2]);
  GetUnitVector(init_univec_normal[crack_count]);
  fclose(fp);
}

void WriteCrackParam(const char *filename, double temp_init_univec_normal[DIMENSION], int crack_count)
{
  FILE *fp;

  fp = fopen(filename, "w");
  assert(fp != NULL);

  fprintf(fp, "%lf\n", nodal_distance[crack_count]);
  fprintf(fp, "%lf %lf %lf\n", temp_init_univec_normal[0], temp_init_univec_normal[1], temp_init_univec_normal[2]);
  fclose(fp);
}

void WriteCrackFrontPoints(const char *filename, int crack_count)
{
  FILE *fp;

  int i;

  fp = fopen(filename, "w");
  assert(fp != NULL);

  fprintf(fp, "%d\n", nnodes[crack_count]);

  for(i = 0; i < nnodes[crack_count]; i++){
    fprintf(fp, "%d %lf %lf %lf\n", i, nodes_coordinate[crack_count][i][0], nodes_coordinate[crack_count][i][1], nodes_coordinate[crack_count][i][2]);
  }
  fclose(fp);
}

void InitializeFlag(int flag[], int total)
{
  int i;
  for(i = 0; i < total; i++){
    flag[i] = 0;
  }
}

/*
 * 入力した２つのベクトルの平均を出力する関数
 */
void GetAverageVector(double a[], double b[], double c[])
{
  int i;
  for(i = 0; i < DIMENSION; i++){
    c[i] = (a[i] + b[i])/2;
  }
}



/*
 * 各前縁点における接線方向ベクトル、法線方向ベクトル、進展方向ベクトルを作成
 * ただし法線方向ベクトルはReadParamで読み取った値を使用
 */
void GetUnivectoratAllPoint(int temp_nnodes, double temp_p_to_p_vector[][DIMENSION], double temp_univec_normal[][DIMENSION], double temp_univec_propa[][DIMENSION], double temp_univec_tangent[][DIMENSION], int crack_count)
{
  int i, j;
  for(j = 0; j < DIMENSION; j++){
    temp_univec_tangent[0][j] = temp_p_to_p_vector[0][j];
    temp_univec_normal[0][j] = init_univec_normal[crack_count][j];
  }
  CrossProduct(temp_univec_normal[0], temp_univec_tangent[0], temp_univec_propa[0]);
  GetUnitVector(temp_univec_propa[0]);


  for(i = 1; i < temp_nnodes - 1; i++){
    for(j = 0; j < DIMENSION; j++){
      temp_univec_normal[i][j] = init_univec_normal[crack_count][j];
    }
    GetAverageVector(temp_p_to_p_vector[i-1], temp_p_to_p_vector[i], temp_univec_tangent[i]);
    CrossProduct(temp_univec_normal[i], temp_univec_tangent[i], temp_univec_propa[i]);
    GetUnitVector(temp_univec_propa[i]);
  }
  for(j = 0; j < DIMENSION; j++){
    temp_univec_tangent[temp_nnodes-1][j] = temp_p_to_p_vector[temp_nnodes-2][j];
    temp_univec_normal[temp_nnodes-1][j] = init_univec_normal[crack_count][j];
  }
  CrossProduct(temp_univec_normal[temp_nnodes-1], temp_univec_tangent[temp_nnodes-1], temp_univec_propa[temp_nnodes-1]);
  GetUnitVector(temp_univec_propa[temp_nnodes-1]);
}

/*
 * き裂前縁点の座標値を用いて始点から終点に向かう方向に１点ずつ点を結ぶベクトルを作成
 */
void GetPointtoPointVector(int temp_nnodes, double temp_nodes_coordinate[][DIMENSION], double temp_p_to_p_vector[][DIMENSION])
{
  int i, j;
  for(i = 0; i < temp_nnodes-1; i++){
    for(j = 0; j < DIMENSION; j++){
      temp_p_to_p_vector[i][j] = temp_nodes_coordinate[i+1][j] - temp_nodes_coordinate[i][j];
      //printf("##temp_p_to_p_vector[%d][%d] = temp_nodes_coordinate[%d][%d] - temp_nodes_coordinate[%d][%d]\n%lf = %lf - %lf\n", i, j, i+1, j, i, j, temp_p_to_p_vector[i][j], temp_nodes_coordinate[i+1][j], temp_nodes_coordinate[i][j]);
    }
    GetUnitVector(temp_p_to_p_vector[i]);
    //printf("##%dto%dvector = (%lf, %lf, %lf)\n", i, i+1, temp_p_to_p_vector[i][0],temp_p_to_p_vector[i][1],temp_p_to_p_vector[i][2]);
  }
}

void GetPairMinDist(int crack_1, int crack_2, int *min_1, int *min_2)
{
  int i,j;
  double min_value;
  double temp_value;
  *min_1 = 0;
  *min_2 = 0;
  min_value = Distance(nodes_coordinate[crack_1][0], nodes_coordinate[crack_2][0]);
  for(i = 0; i < nnodes[crack_1]; i++){
    for(j = 0; j < nnodes[crack_2]; j++){
      temp_value = Distance(nodes_coordinate[crack_1][i], nodes_coordinate[crack_2][j]);
      if(temp_value < min_value){
        min_value = temp_value;
        *min_1 = i;
        *min_2 = j;
      }
    }
  }
  //printf("min_1 = %d, min_2 = %d\n", *min_1, *min_2);
}

#define MINUS 0
#define PLUS 1

int CheckNumberArea(int crack_count, int check_number)
{
  if(nnodes[crack_count] < check_number) return 1;
  if(check_number < 0) return 1;
  return 0;
}

void SplainPoints(double a[], double b[], double c[], int number_of_division)
{
  int i, j;
  double splain_weight = 1.0 / (double)number_of_division;
  double weight;
  for(i = 1; i < number_of_division; i++){
    weight = splain_weight * (double)i;
    printf("weight = %lf\n", weight);
    for(j = 0; j < DIMENSION; j++){
      coalesced_nodes_coordinate[number_of_coalescence_points][j] = a[j] * pow(1.0 - weight, 2.0) + b[j] * 2.0 * weight * (1.0 - weight) + c[j] * pow(weight, 2.0);
    }
    number_of_coalescence_points++;
  }
}

int MakeBridgePoints(int crack_1, int crack_2, int min_1, int min_2, int order_1, int order_2)
{
  int i, j;
  int check_number_area;
  int anchor_point_1;
  double anchor_point_coordinate_1[DIMENSION];
  int splain_anchor_1;
  double splain_anchor_coordinate_1[DIMENSION];
  int anchor_point_2;
  double anchor_point_coordinate_2[DIMENSION];
  int splain_anchor_2;
  double splain_anchor_coordinate_2[DIMENSION];

  double coord_mid_point_of_anchor[DIMENSION];

  double anchor_dist_1;
  double anchor_dist_2;
  double average_nodal_distance;
  int number_of_division_1;
  int number_of_division_2;
  double temp_nodal_distance_1;
  double temp_nodal_distance_2;

  if(order_1 == MINUS){
    anchor_point_1 = min_1 - 5;
    splain_anchor_1 = min_1 - 10;
    check_number_area = CheckNumberArea(crack_1, splain_anchor_1);
    if(check_number_area == 1) return 1;
  }
  if(order_1 == PLUS){
    anchor_point_1 = min_1 + 5;
    splain_anchor_1 = min_1 + 10;
    check_number_area = CheckNumberArea(crack_1, splain_anchor_1);
    if(check_number_area == 1) return 1;
  }
  for(i = 0; i < DIMENSION; i++){
    anchor_point_coordinate_1[i] = nodes_coordinate[crack_1][anchor_point_1][i];
    splain_anchor_coordinate_1[i] = nodes_coordinate[crack_1][splain_anchor_1][i];
  }

  if(order_2 == MINUS){
    anchor_point_2 = min_2 - 5;
    splain_anchor_2 = min_2 - 10;
    check_number_area = CheckNumberArea(crack_2, splain_anchor_2);
    if(check_number_area == 1) return 1;
  }
  if(order_2 == PLUS){
    anchor_point_2 = min_2 + 5;
    splain_anchor_2 = min_2 + 10;
    check_number_area = CheckNumberArea(crack_2, splain_anchor_2);
    if(check_number_area == 1) return 1;
  }
  //printf("order_1 = %d, order_2 = %d\n", order_1, order_2);
  for(i = 0; i < DIMENSION; i++){
    anchor_point_coordinate_2[i] = nodes_coordinate[crack_2][anchor_point_2][i];
    splain_anchor_coordinate_2[i] = nodes_coordinate[crack_2][splain_anchor_2][i];
  }
  for(i = 0; i < DIMENSION; i++){
    coord_mid_point_of_anchor[i] = (nodes_coordinate[crack_1][anchor_point_1][i] + nodes_coordinate[crack_2][anchor_point_2][i])/2;
  }
  anchor_dist_1 = Distance(anchor_point_coordinate_1, coord_mid_point_of_anchor);
  anchor_dist_2 = Distance(anchor_point_coordinate_2, coord_mid_point_of_anchor);
  //printf("anchor_dist_1 = %lf\n", anchor_dist_1);
  //printf("anchor_dist_2 = %lf\n", anchor_dist_2);
  /*
     nodal_distanceはどうするか決めていないのでとりあえず使用した二つのき裂の平均をとる
     そのうちき裂の大きさや曲率に依存して変えるようにするとよさげ？
     VCCMから相互積分法にするに当たって変えなければいけない部分だと思うので
     今は特に考えません
     */
  average_nodal_distance = nodal_distance[crack_1] + nodal_distance[crack_2];
  number_of_division_1 = anchor_dist_1 / average_nodal_distance + 0.5;
  number_of_division_2 = anchor_dist_2 / average_nodal_distance + 0.5;
  temp_nodal_distance_1 = anchor_dist_1 / number_of_division_1;
  temp_nodal_distance_2 = anchor_dist_2 / number_of_division_2;

  number_of_coalescence_points = 0;
  if(order_1 == PLUS){
    for(i = nnodes[crack_1] - 1; i >= splain_anchor_1; i--){
      for(j = 0; j < DIMENSION; j++){
        coalesced_nodes_coordinate[number_of_coalescence_points][j] = nodes_coordinate[crack_1][i][j];
      }
        number_of_coalescence_points++;
    }
  }
  if(order_1 == MINUS){
    for(i = 0; i <= splain_anchor_1; i++){
      for(j = 0; j < DIMENSION; j++){
        coalesced_nodes_coordinate[number_of_coalescence_points][j] = nodes_coordinate[crack_1][i][j];
      }
        number_of_coalescence_points++;
    }
  }
  //printf("number_of_coalescence_points = %d\n", number_of_coalescence_points);
  //printf("number_of_division_1 = %d\n", number_of_division_1);
  //SplainPoints(splain_anchor_coordinate_1, anchor_point_coordinate_1, coord_mid_point_of_anchor, 5 + number_of_division_1);
  //printf("number_of_coalescence_points = %d\n", number_of_coalescence_points);
  for(i = 0; i < DIMENSION; i++){
    coalesced_nodes_coordinate[number_of_coalescence_points][i] = coord_mid_point_of_anchor[i];
  }
  number_of_coalescence_points++;
  //printf("number_of_coalescence_points = %d\n", number_of_coalescence_points);
  //printf("number_of_division_2 = %d\n", number_of_division_2);
  SplainPoints(coord_mid_point_of_anchor , anchor_point_coordinate_2, splain_anchor_coordinate_2, 100 + number_of_division_2);
  //printf("number_of_coalescence_points = %d\n", number_of_coalescence_points);
  if(order_2 == PLUS){
    for(i = splain_anchor_2; i < nnodes[crack_2]; i++){
      for(j = 0; j < DIMENSION; j++){
        coalesced_nodes_coordinate[number_of_coalescence_points][j] = nodes_coordinate[crack_2][i][j];
      }
        number_of_coalescence_points++;
    }
  }
  if(order_2 == MINUS){
    for(i = splain_anchor_2; i >= 0; i--){
      for(j = 0; j < DIMENSION; j++){
        coalesced_nodes_coordinate[number_of_coalescence_points][j] = nodes_coordinate[crack_2][i][j];
      }
        number_of_coalescence_points++;
    }
  }
  printf("number_of_coalescence_points = %d\n", number_of_coalescence_points);
  return 0;
}

void TempResisterCoalescedPoints(int number_of_coalescence)
{
  int i, j;
  temp_number_of_coalescence_points[number_of_coalescence] = number_of_coalescence_points;
  for(i = 0; i < temp_number_of_coalescence_points[number_of_coalescence]; i++){
    for(j = 0; j < DIMENSION; j++){
      temp_coalescence_points_coordinte[number_of_coalescence][i][j] = coalesced_nodes_coordinate[i][j];
    }
  }
}

/*
   p_tp_p_vectorのベクトル番号（左から二番目の配列）は
   その番号の点から次の点へのベクトルを示している
   */
int GetBridge(int crack_1, int crack_2, int min_1, int min_2)
{
  int i;
  int coalescence_flag = 0;
  int number_of_coalescence = 0;
  double temp_normal_crack_1[DIMENSION];
  double temp_normal_crack_2[DIMENSION];
  double temp_inner_product;
  double minus_p_to_p_vector[2][DIMENSION];
  //printf("p_tp_p_vector[%d][%d] = %lf %lf %lf\np_tp_p_vector[%d][%d] = %lf %lf %lf\n", crack_1, min_1, p_to_p_vector[crack_1][min_1][0],p_to_p_vector[crack_1][min_1][1], p_to_p_vector[crack_1][min_1][2], crack_2, min_2, p_to_p_vector[crack_2][min_2][0], p_to_p_vector[crack_2][min_2][1], p_to_p_vector[crack_2][min_2][2]);

  if(min_1 != nnodes[crack_1]-1 && min_2 != nnodes[crack_2]-1){
  CrossProduct(univec_propa[crack_1][min_1], p_to_p_vector[crack_1][min_1], temp_normal_crack_1);
  CrossProduct(univec_propa[crack_2][min_2], p_to_p_vector[crack_2][min_2], temp_normal_crack_2);
  temp_inner_product = InnerProduct(temp_normal_crack_1, temp_normal_crack_2);
  printf("temp_inner_product = %lf\n", temp_inner_product);
  if(temp_inner_product < 0) coalescence_flag = MakeBridgePoints(crack_1, crack_2, min_1, min_2, PLUS, PLUS);
  if(coalescence_flag == 0){
    printf("coalescence_flag == 1\n");
    TempResisterCoalescedPoints(number_of_coalescence);
    number_of_coalescence++;
  }
    coalescence_flag = 0;
  }

  for(i = 0; i < DIMENSION; i++){
    minus_p_to_p_vector[0][i] = -p_to_p_vector[crack_1][min_1-1][i];
    minus_p_to_p_vector[1][i] = -p_to_p_vector[crack_2][min_2-1][i];
  }

  if(min_1 != 0 && min_2 != nnodes[crack_2]-1){
  CrossProduct(univec_propa[crack_1][min_1], minus_p_to_p_vector[0], temp_normal_crack_1);
  CrossProduct(univec_propa[crack_2][min_2], p_to_p_vector[crack_2][min_2], temp_normal_crack_2);
  temp_inner_product = InnerProduct(temp_normal_crack_1, temp_normal_crack_2);
  printf("temp_inner_product = %lf\n", temp_inner_product);
  if(temp_inner_product < 0) coalescence_flag = MakeBridgePoints(crack_1, crack_2, min_1, min_2, MINUS, PLUS);
  if(coalescence_flag == 0){
    printf("coalescence_flag == 1\n");
    TempResisterCoalescedPoints(number_of_coalescence);
    number_of_coalescence++;
  }
    coalescence_flag = 0;
  }

  if(min_1 != nnodes[crack_1]-1 && min_2 != 0){
  CrossProduct(univec_propa[crack_1][min_1], p_to_p_vector[crack_1][min_1], temp_normal_crack_1);
  CrossProduct(univec_propa[crack_2][min_2], minus_p_to_p_vector[1], temp_normal_crack_2);
  temp_inner_product = InnerProduct(temp_normal_crack_1, temp_normal_crack_2);
  printf("temp_inner_product = %lf\n", temp_inner_product);
  if(temp_inner_product < 0) coalescence_flag = MakeBridgePoints(crack_1, crack_2, min_1, min_2, PLUS, MINUS);
  if(coalescence_flag == 0){
    printf("coalescence_flag == 1\n");
    TempResisterCoalescedPoints(number_of_coalescence);
    number_of_coalescence++;
  }
    coalescence_flag = 0;
  }


  if(min_1 != 0 && min_2 != 0){
  CrossProduct(univec_propa[crack_1][min_1], minus_p_to_p_vector[0], temp_normal_crack_1);
  CrossProduct(univec_propa[crack_2][min_2], minus_p_to_p_vector[1], temp_normal_crack_2);
  temp_inner_product = InnerProduct(temp_normal_crack_1, temp_normal_crack_2);
  printf("temp_inner_product = %lf\n", temp_inner_product);
  if(temp_inner_product < 0) coalescence_flag = MakeBridgePoints(crack_1, crack_2, min_1, min_2, MINUS, MINUS);
  if(coalescence_flag == 0){
    printf("coalescence_flag == 1\n");
    TempResisterCoalescedPoints(number_of_coalescence);
    number_of_coalescence++;
  }
    coalescence_flag = 0;
  }
  return number_of_coalescence;
}

void UpdateCrackNumbers(int crack_1, int crack_2, int number_of_coalescence)
{
  int i, j;
  int crack_count;
  int temp_coalescence_crack[NUMBER_OF_CRACKS][2];
  int temp_coalescence_crack_flag[2][NUMBER_OF_CRACKS];
  for(i = 0; i < 2; i++){
    InitializeFlag(temp_coalescence_crack_flag[i], NUMBER_OF_CRACKS);
  }

  assert(number_of_coalescence != 0);
  /*
     二つのき裂が表面付近で合体し、一つのき裂になる
     */
  if(number_of_coalescence == 1){
    nnodes[crack_1] = temp_number_of_coalescence_points[0];
    for(i = 0; i < temp_number_of_coalescence_points[0]; i++){
      for(j = 0; j < DIMENSION; j++){
        nodes_coordinate[crack_1][i][j] = temp_coalescence_points_coordinte[0][i][j];
      }
    }
    for(j = 0; j < coalescence_count; j++){
      if(coalescence_crack[j][0] == crack_2){
        coalescence_crack[j][0] = crack_1;
      }
      if(coalescence_crack[j][1] == crack_2){
        coalescence_crack[j][1] = crack_1;
      }
    }
    for(crack_count = crack_2; crack_count < cracks - 1; crack_count++){
      for(i = 0; i < nnodes[crack_count+1]; i++){
        for(j = 0; j < DIMENSION; j++){
          nodes_coordinate[crack_count][i][j] = nodes_coordinate[crack_count+1][i][j];
          univec_normal[crack_count][i][j] = univec_normal[crack_count+1][i][j];
          univec_propa[crack_count][i][j] = univec_propa[crack_count+1][i][j];
          univec_tangent[crack_count][i][j] = univec_tangent[crack_count+1][i][j];
        }
      }
      for(i = 0; i < nnodes[crack_count+1] - 1; i++){
        for(j = 0; j < DIMENSION; j++){
          p_to_p_vector[crack_count][i][j] = p_to_p_vector[crack_count+1][i][j];
        }
      }
      nnodes[crack_count] = nnodes[crack_count+1];
      nodal_distance[crack_count] = nodal_distance[crack_count+1];
      for(j = 0; j < coalescence_count; j++){
        if(coalescence_crack[j][0] == crack_count){
          temp_coalescence_crack[j][0] = coalescence_crack[j][0] - 1;
          temp_coalescence_crack_flag[0][j] = 1;
        }
        if(coalescence_crack[j][1] == crack_count){
          temp_coalescence_crack[j][1] = coalescence_crack[j][1] - 1;
          temp_coalescence_crack_flag[1][j] = 1;
        }
      }
    }
    for(i = 0; i < coalescence_count; i++){
      if(temp_coalescence_crack_flag[0][i] == 1){
        coalescence_crack[i][0] = temp_coalescence_crack[i][0];
      }
      if(temp_coalescence_crack_flag[1][i] == 1){
        coalescence_crack[i][1] = temp_coalescence_crack[i][1];
      }
    }
    cracks--;
  }

  /*
     合体の様式として二つのき裂が中央部付近で合体して
     二つのき裂になるための
     クライテリアを用意してないのでまだ作り込みません*/
  /*
     if(number_of_coalescence == 2){
     nnodes[crack_1] = temp_number_of_coalescence_points[0];
     for(i = 0; i < temp_number_of_coalescence_points[0]; i++){
     for(j = 0; j < DIMENSION; j++){
     nodes_coordinate[crack_1][i][j] = temp_coalescence_points_coordinte[0][i][j];
     }
     }
     }*/
}

void BridgePoints(int crack_1, int crack_2)
{
  int *ptr_min_1, *ptr_min_2;
  int min_1, min_2;
  int number_of_coalescence;
  ptr_min_1 = &min_1;
  ptr_min_2 = &min_2;
  GetPairMinDist(crack_1, crack_2, ptr_min_1, ptr_min_2);
  printf("min_1 = %d, min_2 = %d\n", min_1, min_2);
  number_of_coalescence = GetBridge(crack_1, crack_2, min_1, min_2);
  printf("number_of_coalescence = %d\n", number_of_coalescence);
  UpdateCrackNumbers(crack_1, crack_2, number_of_coalescence);
}


/*自己き裂合体はまだ未実装
  void BridgeSelfPoints(int crack)
  {

  }
  */

void PerformCommand()
{
  int i;
  int temp_crack_num[2];
  for(i = 0; i < coalescence_count; i++){
    temp_crack_num[0] = coalescence_crack[i][0];
    temp_crack_num[1] = coalescence_crack[i][1];
    if(coalescence_crack[i][0] != coalescence_crack[i][1]){
      BridgePoints(coalescence_crack[i][0], coalescence_crack[i][1]);
    } /*else {自己き裂合体はまだ未実装
        BridgeSelfPoints(coalescence_crack[i][0]);
        }*/
  }
}

int main(int argc, char *argv[])
{
  int i, j;
  ReadCoalescenceFlag(argv[1]);
  int cracks;
  cracks = atoi(argv[2]);
  int crack_count;
  char front_points_filename[30];
  char param_filename[30];
  double temp_init_univec_normal[DIMENSION];

  for(crack_count = 0; crack_count < cracks; crack_count++){
    sprintf(front_points_filename, "%s%d", argv[3], crack_count);
    ReadCrackFrontPoints(front_points_filename, crack_count);
    sprintf(param_filename, "%s%d", argv[4], crack_count);
    ReadCrackParam(param_filename, crack_count);
    GetPointtoPointVector(nnodes[crack_count], nodes_coordinate[crack_count], p_to_p_vector[crack_count]);
    GetUnivectoratAllPoint(nnodes[crack_count], p_to_p_vector[crack_count], univec_normal[crack_count], univec_propa[crack_count], univec_tangent[crack_count], crack_count);
  }


  PerformCommand();


  for(crack_count = 0; crack_count < cracks; crack_count++){
    GetPointtoPointVector(nnodes[crack_count], nodes_coordinate[crack_count], p_to_p_vector[crack_count]);
    GetUnivectoratAllPoint(nnodes[crack_count], p_to_p_vector[crack_count], univec_normal[crack_count], univec_propa[crack_count], univec_tangent[crack_count], crack_count);
    for(i = 0; i < nnodes[crack_count]; i++){
      for(j = 0; j < DIMENSION; j++){
        temp_init_univec_normal[j] += univec_normal[crack_count][i][j];
      }
    }
    for(j = 0; j < DIMENSION; j++){
      temp_init_univec_normal[j] /= nnodes[crack_count];
    }
    sprintf(param_filename, "%s%d_refined", argv[4], crack_count);
    WriteCrackParam(param_filename, temp_init_univec_normal, crack_count);
    sprintf(front_points_filename, "%s%d_refined", argv[3], crack_count);
    WriteCrackFrontPoints(front_points_filename, crack_count);
  }

  return 0;
}
