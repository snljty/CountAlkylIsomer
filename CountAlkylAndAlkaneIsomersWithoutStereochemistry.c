# include <stdio.h>

/*****************************************************************************
 *  Using the generating function, Burnside's lemma and Polya's enumeration  *
 *  theory to compute the amount of isomers for alkyl and alkane without     *
 *  considering stereochemistry.                                             *
 *                                                                           *
 *  The generating function for alkyl is:                                    *
 *  A(x) = 1 + A_1 x + A_2 x^2 + A_3 x^3 + ... ,                             *
 *  in which A_m is the amount of alkyl isomers with m carbons.              *
 *  Based on Polya's theory,                                                 *
 *  A(x) = 1 + x * [A(x)^3 + 3 A(x) A(x^2) + 2 A(x^3)] / 6 .                 *
 *  Using A_1, A_2, ..., A_m, we can figure out A_(m + 1) .                  *
 *                                                                           *
 *  For alkane with one indexed carbon atom,                                 *
 *  P(x) = P_1 x + P_2 x^2 + P_3 x^3 + ... ,                                 *
 *  and we have:                                                             *
 *  P(x) = x * [A(x)^4 + 6*A(x)^2*A(x^2) + 3*A(x^2)^2 + 8*A(x)*A(x^3)        *
 *              + 6A(x^4)] / 24 .                                            *
 *  For alkane with one indexed carbon-carbon bond,                          *
 *  Q(x) = Q_1 x + Q_2 x^2 + Q_3 x^3 + ... ,                                 *
 *  and we have:                                                             *
 *  Q(x) = [A(x)^2 + A(x^2)] / 2 - A(x) .                                    *
 *  Finally, for alkane,                                                     *
 *  C(x) = C_1 x + C_2 x^2 + C_3 x^3 + ... ,                                 *
 *  and we also have:                                                        *
 *  C(x) = P(x) - Q(x) + A(x^2) - 1 .                                        *
 *****************************************************************************/

# define MAX_CARBON 30
# define CALC_CARBON 20
int A[MAX_CARBON + 1] = {0};
int P[MAX_CARBON + 1] = {0}, Q[MAX_CARBON + 1] = {0};
int C[MAX_CARBON + 1] = {0};

void TryCalc_A_m(int m);

void Calc_A_m(int m);

void TryCalc_P_m(int m);

void Calc_P_m(int m);

void TryCalc_Q_m(int m);

void Calc_Q_m(int m);

void TryCalc_C_m(int m);

void Calc_C_m(int m);


int main()
{
  int num = 0;

  Calc_A_m(0);
  for (num = 1; num <= CALC_CARBON; ++ num)
  {
    Calc_A_m(num);
    Calc_P_m(num);
    Calc_Q_m(num);
    Calc_C_m(num);
  }
  num = 1;
  printf("CH%d-: %d\n", 2 * num + 1, A[num]);
  for (num = 2; num <= CALC_CARBON; ++ num)
    printf("C%dH%d-: %d\n", num, 2 * num + 1, A[num]);
  puts("");
  num = 1;
  printf("CH%d: %d\n", 2 * num + 2, A[num]);
  for (num = 2; num <= CALC_CARBON; ++ num)
    printf("C%dH%d: %d\n", num, 2 * num + 2, C[num]);

  return 0;
}


void TryCalc_A_m(int m)
{
  /*  m >= 0 for A[m]  */
  int i = 0;

  for (i = 0; i <= m; ++ i)
    if (! A[i])
      Calc_A_m(i);
  return;
}

void Calc_A_m(int m)
{
  /*  m >= 0 for A[m]  */
  int t = 0, ans = 0;
  int i = 0, j = 0;

  /*  A[0] = 1  */
  if (! m)
  {
    A[m] = 1;
    return;
  }
  /*  term 1: x * A(x)^3  */
  t = 0;
  for (i = 0; i < m; ++ i)
    for (j = 0; i + j < m; ++ j)
      t += A[i] * A[j] * A[m - 1 - i - j];
  ans += t;
  /*  term 2: x * 3 * A(x) * A(x^2)  */
  t = 0;
  for (j = 0; 2 * j < m; ++ j)
    t += A[m - 1 - 2 * j] * A[j];
  ans += 3 * t;
  /* term 3: x * 2 * A(x^3)  */
  t = 0;
  if (! ((m - 1) % 3))
    t = A[(m - 1) / 3];
  ans += 2 * t;
  /*  divided by 6  */
  ans /= 6;
  A[m] = ans;
  return;
}

void TryCalc_P_m(int m)
{
  /*  m >= 1 for P[m]  */
  int i = 0;

  for (i = 0; i <= m; ++ i)
    if (! A[i])
      Calc_A_m(i);
  if (! P[m])
    Calc_P_m(m);
  return;
}

void Calc_P_m(int m)
{
  /*  m >= 1 for P[m]  */
  int t = 0, ans = 0;
  int i = 0, j = 0, k = 0;

  /*  term 1: x * A(x)^4  */
  t = 0;
  for (i = 0; i < m; ++ i)
    for (j = 0; i + j < m; ++ j)
      for (k = 0; i + j + k < m; ++ k)
        t += A[i] * A[j] * A[k] * A[m - 1 - i - j - k];
  ans += t;
  /*  term 2: x * 6 * A(x)^2 * A(x^2)  */
  t = 0;
  for (k = 0; 2 * k < m; ++ k)
    for (i = 0; i + 2 * k < m; ++ i)
      t += A[i] * A[m - 1 - i - 2 * k] * A[k];
  ans += 6 * t;
  /*  term 3: x * 3 * A(x^2)^2  */
  t = 0;
  if (! ((m - 1) % 2))
    for (i = 0; 2 * i < m; ++ i)
      t += A[i] * A[(m - 1 - 2 * i) / 2];
  ans += 3 * t;
  /*  term 4: x * 8 * A(x) * A(x^3)  */
  t = 0;
  for (j = 0; 3 * j < m; ++ j)
    t += A[m - 1 - 3 * j] * A[j];
  ans += 8 * t;
  /*  term 5: x * 6 * A(x^4)  */
  t = 0;
  if (! ((m - 1) % 4))
    t = A[(m - 1) / 4];
  ans += 6 * t;
  /*  divided by 24  */
  ans /= 24;
  P[m] = ans;
  return;
}

void TryCalc_Q_m(int m)
{
  /*  m >= 1 for Q[m]  */
  int i = 0;

  for (i = 0; i <= m; ++ i)
    if (! A[i])
      Calc_A_m(i);
  if (! Q[m] && m != 1)
    Calc_Q_m(m);
  return;
}

void Calc_Q_m(int m)
{
  /*  m >= 1 for Q[m]  */
  int t = 0, ans = 0;
  int i = 0;

  if (m == 1)
  {
    Q[m] == 0;
    return;
  }
  /*  term 1: A(x)^2  */
  t = 0;
  for (i = 0; i <= m; ++ i)
    t += A[i] * A[m - i];
  ans += t;
  /*  term 2: A(x^2)  */
  t = 0;
  if (! (m % 2))
    t = A[m / 2];
  ans += t;
  /*  divided by 2  */
  ans /= 2;
  /*  term 3: - A(x)  */
  ans -= A[m];
  Q[m] = ans;
  return;
}

void TryCalc_C_m(int m)
{
  /*  m >= 1 for C[m]  */
  int i = 1;

  if (! A[0])
    Calc_A_m(0);
  if (! A[1])
    Calc_A_m(1);
  if (! P[1])
    Calc_P_m(1);
  for (i = 2; i <= m; ++ i)
  {
    if (! A[i])
      Calc_A_m(i);
    if (! P[i])
      Calc_P_m(i);
    if (! Q[i])
      Calc_Q_m(i);
  }
  return;
}

void Calc_C_m(int m)
{
  /*  m >= 1 for C[m]  */
  int ans = 0;
  if (! (m % 2))
    ans = P[m] - Q[m] + A[m / 2];
  else
    ans = P[m] - Q[m];
  C[m] = ans;
  return;
}

