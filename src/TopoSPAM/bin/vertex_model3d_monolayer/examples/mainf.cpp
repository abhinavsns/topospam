#include <iostream>

#include "Bond-inline.h"

int main()
{
  Vector3d v1(0,0,1);
  Vector3d v2(1,0,0);
  Vector3d v3(0,2,1);

  Matrix3d Q3d;
  Vector3d n_t;
  double area;

  std::tie(Q3d, n_t, area) = calculate_triangle_elongation_tensor(v1, v2, v3);

  std::cout << Q3d << std::endl;
  std::cout << n_t << std::endl;
  std::cout << area << std::endl;

  auto [val, vec] = calculate_3x3_tensor_eigenvalues_and_eigenvector(Q3d, 1);

  std::cout << "val:\n" << val << std::endl;
  std::cout << "vec:\n" << vec << std::endl;

  Vector3d nn = (v2-v1).cross(v3-v1);
  auto rot_Q3d = rotate_tensors_in_plane(Q3d, nn/nn.norm());

  std::cout << "rot_Q3d:\n" << rot_Q3d << std::endl;

  auto [val2, vec2] = calculate_3x3_tensor_eigenvalues_and_eigenvector(rot_Q3d, 1);

  std::cout <<  "val2:\n" << val2 << std::endl;
  std::cout <<  "vec2:\n" << vec2 << std::endl;
  return 0;
}
