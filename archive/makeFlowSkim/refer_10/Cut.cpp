#include "Cut.h"

Cut::Cut() : m_enabled(false), m_value(0) {
  // skip
}

bool Cut::enabled() {
  return m_enabled;
}

void Cut::set(double value) {
  // set the cut
  m_value = value;

  // enable the cut
  m_enabled = true;
}

bool Cut::operator>(const double& rhs) {
  return m_value > rhs;
}

bool Cut::operator<(const double &rhs)
{
  return m_value < rhs;
}

bool Cut::operator>=(const double &rhs)
{
  return m_value >= rhs;
}

bool Cut::operator<=(const double &rhs)
{
  return m_value <= rhs;
}

bool Cut::operator==(const double &rhs)
{
  return m_value == rhs;
}

bool Cut::operator!=(const double &rhs)
{
  return m_value != rhs;
}