module global
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Using the generating function, Burnside's lemma and Polya's enumeration  !
  !  theory to compute the amount of isomers for alkyl and alkane without     !
  !  considering stereochemistry.                                             !
  !                                                                           !
  !  The generating function for alkyl is:                                    !
  !  A(x) = 1 + A_1 x + A_2 x^2 + A_3 x^3 + ... ,                             !
  !  in which A_m is the amount of alkyl isomers with m carbons.              !
  !  Based on Polya's theory,                                                 !
  !  A(x) = 1 + x * [A(x)^3 + 3 A(x) A(x^2) + 2 A(x^3)] / 6 .                 !
  !  Using A_1, A_2, ..., A_m, we can figure out A_(m + 1) .                  !
  !                                                                           !
  !  For alkane with one indexed carbon atom,                                 !
  !  P(x) = P_1 x + P_2 x^2 + P_3 x^3 + ... ,                                 !
  !  and we have:                                                             !
  !  P(x) = x * [A(x)^4 + 6*A(x)^2*A(x^2) + 3*A(x^2)^2 + 8*A(x)*A(x^3)        !
  !              + 6A(x^4)] / 24 .                                            !
  !  For alkane with one indexed carbon-carbon bond,                          !
  !  Q(x) = Q_1 x + Q_2 x^2 + Q_3 x^3 + ... ,                                 !
  !  and we have:                                                             !
  !  Q(x) = [A(x)^2 + A(x^2)] / 2 - A(x) .                                    !
  !  Finally, for alkane,                                                     !
  !  C(x) = C_1 x + C_2 x^2 + C_3 x^3 + ... ,                                 !
  !  and we also have:                                                        !
  !  C(x) = P(x) - Q(x) + A(x^2) - 1 .                                        !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  implicit none
  integer, PARAMETER :: MAX_CARBON = 30
  integer, PARAMETER :: CALC_CARBON = 20
  integer :: A(0: MAX_CARBON)
  integer :: P(0: MAX_CARBON), Q(0: MAX_CARBON)
  integer :: C(0: MAX_CARBON)
end module

program main
  use global
  implicit none

  integer :: NUM = 0

  do NUM = 0, MAX_CARBON, 1
    A(NUM) = 0
    P(NUM) = 0
    Q(NUM) = 0
    C(NUM) = 0
  end do

  call Calc_A_m(0)
  do NUM = 1, CALC_CARBON, 1
    call Calc_A_m(NUM)
    call Calc_P_m(NUM)
    call Calc_Q_m(NUM)
    call Calc_C_m(NUM)
  end do
  write(*, "(A, A, A, I2, A, I8)") "C", "  ", "H", (2 * 1 + 1), " -: ", A(1)
  do NUM = 2, CALC_CARBON, 1
    write(*, "(A, I2, A, I2, A, I8)") "C", NUM, "H", (2 * NUM + 1), " -: ", A(NUM)
  end do
  write(*, *) ""
  write(*, "(A, A, A, I2, A, I8)") "C", "  ", "H", (2 * 1 + 2), " : ", C(1)
  do NUM = 2, CALC_CARBON, 1
    write(*, "(A, I2, A, I2, A, I8)") "C", NUM, "H", (2 * NUM + 2), " : ", C(NUM)
  end do

  stop
end program

subroutine TryCalc_A_m(M)
  use global
  implicit none

  integer :: M
  ! M >= 0 for A(M)
  integer :: I = 0

  do I = 0, m, 1
    if(A(I) == 0) then
      call Calc_A_m(i)
    end if
  end do
  return
end subroutine

subroutine Calc_A_m(M)
  use global
  implicit none

  integer :: M
  !  M >= 0 for A(M)
  integer :: T = 0, ANS = 0
  integer :: I = 0, J = 0

  !  A(0) = 1
  ANS = 0
  if(M == 0) then
    A(M) = 1
    return
  end if
  !  term 1: x * A(x)^3
  T = 0
  I = 0
  do while(I < M)
    J = 0
    do while(I + J < M)
      T = T + A(I) * A(J) * A(M - 1 - I - J)
      J = J + 1
    end do
    I = I + 1
  end do
  ANS = ANS + T
  !  term 2: x * 3 * A(x) * A(x^2)
  T = 0
  J = 0
  do while(2 * J < M)
    T = T + A(M - 1 - 2 * J) * A(J)
    J = J + 1
  end do
  ANS = ANS + 3 * T
  ! term 3: x * 2 * A(x^3)
  T = 0
  if(mod((M - 1), 3) == 0) then
    T = A((M - 1) / 3)
  end if
  ANS = ANS + 2 * T
  !  divided by 6
  ANS = ANS / 6
  A(M) = ANS
  return
end subroutine

subroutine TryCalc_P_m(M)
  use global
  implicit none

  integer :: M
  !  M >= 1 for P(M)
  integer :: I = 0

  do I = 0, M, 1
    if(A(I) == 0) then
      call Calc_A_m(I)
    end if
  end do
  if(P(M) == 0) then
    call Calc_P_m(M)
  end if
  return
end subroutine

subroutine Calc_P_m(M)
  use global
  implicit none

  integer :: M

  !  M >= 1 for P(M)
  integer :: T = 0, ANS = 0
  integer :: I = 0, J = 0, K = 0

  ANS = 0
  !  term 1: x * A(x)^4
  T = 0
  I = 0
  do while(I < M)
    J = 0
    do while(I + J < M)
      K = 0
      do while(I + J + K < M)
        T = T + A(I) * A(J) * A(K) * A(M - 1 - I - J - K)
        K = K + 1
      end do
      J = J + 1
    end do
    I = I + 1
  end do
  ANS = ANS + T
  !  term 2: x * 6 * A(x)^2 * A(x^2)
  T = 0
  K = 0
  do while(2 * K < M)
    I = 0
    do while(I + 2 * K < M)
      T = T + A(I) * A(M - 1 - I - 2 * K) * A(K)
      I = I + 1
    end do
    K = K + 1
  end do
  ANS = ANS + 6 * T
  !  term 3: x * 3 * A(x^2)^2
  T = 0
  if(mod((M - 1), 2) == 0) then
    I = 0
    do while(2 * I < M)
      T = T + A(I) * A((M - 1) / 2 - I)
      I = I + 1
    end do
  end if
  ANS = ANS + 3 * T
  !  term 4: x * 8 * A(x) * A(x^3)
  T = 0
  J = 0
  do while(3 * J < M)
    T = T + A(M - 1 - 3 * J) * A(J)
    J = J + 1
  end do
  ANS = ANS + 8 * T
  !  term 5: x * 6 * A(x^4)
  T = 0
  if(mod((M - 1), 4) == 0) then
    T = A((M - 1) / 4)
  end if
  ANS = ANS + 6 * T
  !  divided by 24
  ANS = ANS / 24
  P(M) = ANS
  return
end subroutine

subroutine TryCalc_Q_m(M)
  use global
  implicit none

  integer :: M

  !  M >= 1 for Q(M)
  integer :: I = 0

  do I = 0, M, 1
    if(A(I) == 0) then
      call Calc_A_m(I)
    end if
  end do
  if(Q(M) == 0 .AND. M /= 1) then
    call Calc_Q_m(M)
  end if
  return
end subroutine

subroutine Calc_Q_m(M)
  use global
  implicit none

  integer :: M

  !  M >= 1 for Q(M)
  integer :: T = 0, ANS = 0
  integer :: I = 0

  if(M == 1) then
    Q(M) = 0
    return
  end if
  ANS = 0
  !  term 1: A(x)^2
  T = 0
  do I = 0, M, 1
    T = T + A(I) * A(M - I)
  end do
  ANS = ANS + T
  !  term 2: A(x^2)
  T = 0
  if(mod(M, 2) == 0) then
    T = A(M / 2)
  end if
  ANS = ANS + T
  !  divided by 2
  ANS = ANS / 2
  !  term 3: - A(x)
  ANS = ANS - A(M)
  Q(M) = ANS
  return
end subroutine

subroutine TryCalc_C_m(M)
  use global
  implicit none

  integer :: M

  !  M >= 1 for C(M)
  integer :: I = 1

  if(A(0) == 0) then
    call Calc_A_m(0)
  end if
  if(A(1) == 0) then
    call Calc_A_m(1)
  end if
  if(P(1) == 0) then
    call Calc_P_m(1)
  end if
  do I = 2, M, 1
    if(A(I) == 0) then
      call Calc_A_m(I)
    end if
    if(P(I) == 0) then
      call Calc_P_m(I)
    end if
    if(Q(I) == 0) then
      call Calc_Q_m(I)
    end if
  end do
  return
end subroutine

subroutine Calc_C_m(M)
  use global
  implicit none

  integer :: M

  !  M >= 1 for C(M)
  integer :: ANS = 0
  ANS = 0
  if(mod(M, 2) == 0) then
    ANS = P(M) - Q(M) + A(M / 2)
  ELSE
    ANS = P(M) - Q(M)
  end if
  C(M) = ANS
  return
end subroutine

