
int nsamples = 10;
void gausTest()
{
  double mean = double(nsamples / 2);
  double sigma = double(nsamples);
  double w0 = TMath::Gaus(mean, mean, sigma, false);
  printf("W0 %.0f mean %.0f sigma %.f \n", w0, mean, sigma);

  for (int i = 0; i <= nsamples; ++i)
  {
    int ifreq = i - nsamples / 2;
    double x = double(ifreq);
    double w = TMath::Gaus(x, mean, sigma, false);
    printf("  i %i x %.0f w %E \n", i, x, w);
  }

  // test complex

  const double PI = std::acos(-1);             // or std::numbers::pi in C++20
  std::complex<double> z5 = std::exp(1i * PI); // Euler's formula
  std::cout << "exp(i * pi) = " << z5 << '\n';

  std::complex<double> z1 = 3. + 1i * PI / 2.;

  printf(" real %f imag %f mag %f phase %f \n", z1.real(), z1.imag(), std::abs(z1), std::arg(z1));

  std::complex<double> z2 = 1. / z1;
  printf(" real %f imag %f mag %f phase %f \n", z2.real(), z2.imag(), std::abs(z2), std::arg(z2));

  std::complex<double> z3 = std::conj(z1) / (std::norm(z1));

  printf(" real %f imag %f mag %f phase %f \n", z3.real(), z3.imag(), std::abs(z3), std::arg(z3));

  double phase = 10;
  std::complex<double> shift = std::exp(phase * 1i);
  std::cout << " phase " << phase << "  " << shift << '\n';
}
