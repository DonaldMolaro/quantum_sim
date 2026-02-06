#pragma once

struct RegisterLayout {
  int target_start;
  int target_end;
  int control_start;
  int control_end;
};

inline RegisterLayout make_shor_layout(int n_t, int n_c)
{
  // Target register is lower bits; control register follows.
  RegisterLayout layout;
  layout.target_start = 0;
  layout.target_end = n_t - 1;
  layout.control_start = n_t;
  layout.control_end = n_t + n_c - 1;
  return layout;
}
