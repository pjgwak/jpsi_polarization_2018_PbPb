#pragma once

class Cut {
  public:
    /**
     * @brief Construct a new uninitialized cut
     */
    Cut();

    /**
     * @brief Set the cut value
     */
    void set(double value);

    /**
     * @brief Return trun if value of this cut was set
     */
    bool enabled();

    /**
     * @brief Overloaded comparison operatiors
     */
    bool operator>(const double &rhs);
    bool operator<(const double &rhs);
    bool operator>=(const double &rhs);
    bool operator<=(const double &rhs);
    bool operator==(const double &rhs);
    bool operator!=(const double &rhs);

  protected:
    /**
     * @brief Indicates whether cut was initialized
     */
    bool m_enabled;

    /**
     * @brief The cut vale
     */
    double m_value;
};