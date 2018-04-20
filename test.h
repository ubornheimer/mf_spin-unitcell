#ifndef TEST_H
#define TEST_H


class test
{
    public:
        /** Default constructor */
        test();
        /** Default destructor */
        virtual ~test();
        /** Access m_Counter
         * \return The current value of m_Counter
         */
        unsigned int GetCounter() { return m_Counter; }
        /** Set m_Counter
         * \param val New value to set
         */
        void SetCounter(unsigned int val) { m_Counter = val; }
    protected:
    private:
        unsigned int m_Counter; //!< Member variable "m_Counter"
};

#endif // TEST_H
