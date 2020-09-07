#include <cppunit/extensions/HelperMacros.h>
#include "../utils.hpp"
class TestUtils : public CppUnit::TestFixture {
    private:

    public:
      static CppUnit::Test *suite()
      {
        CppUnit::TestSuite *suiteOfTests = new CppUnit::TestSuite( "ComplexNumberTest" );
        suiteOfTests->addTest( new CppUnit::TestCaller<TestUtils>( 
                                      "testCircleSectionParameters", 
                                      &TestUtils::testCircleSectionParameters ) );
        suiteOfTests->addTest( new CppUnit::TestCaller<TestUtils>( 
                                      "testRectSectionParameters", 
                                      &TestUtils::testRectSectionParameters ) );
        suiteOfTests->addTest( new CppUnit::TestCaller<TestUtils>( 
                                      "testTrapSectionParameters", 
                                      &TestUtils::testTrapSectionParameters ) );
        return suiteOfTests;
      }
      void setUp()
      {

      }

      void tearDown() 
      {
      
      }
      void testCircleSectionParameters()
      {
        ChannelShape channel;
        channel.CrossSectionShape = CIRC;
        channel.Diameter = 10.000;
        double YY = 5.000;
        double A = 0;
        double P = 0;
        double TW = 0;
        //Half empty
        CalculateSectionFlowParameters(channel,YY,&A,&P,&TW);
        CPPUNIT_ASSERT( A > 39.1 &&  A < 39.3); // Half empty should be A=39.275
        CPPUNIT_ASSERT( P > 15.6 &&  P < 15.8); // Half empty should be P=15.71
        CPPUNIT_ASSERT( TW > 9.9 &&  TW < 10.1); // Half empty should be TW=10.000
        //Full
        YY = 10.000;
        CalculateSectionFlowParameters(channel,YY,&A,&P,&TW);
        CPPUNIT_ASSERT( A > 78.4 &&  A < 78.6); // Full should be A=78.55
        CPPUNIT_ASSERT( P > 31.3 &&  P < 31.5 ); // Full should be P=31.42
        CPPUNIT_ASSERT( TW > 0 &&  TW < 0.02); // Full should be TW=0.000
      }
      void testRectSectionParameters()
      {
        ChannelShape channel;
        channel.CrossSectionShape = RECT;
        channel.Width = 10.000;
        double YY = 5.000;
        double A = 0;
        double P = 0;
        double TW = 0;
        //5 deep
        CalculateSectionFlowParameters(channel,YY,&A,&P,&TW);
        CPPUNIT_ASSERT( A > 49.9 &&  A < 50.1); // should be A=50.00
        CPPUNIT_ASSERT( P > 19.9 &&  P < 20.1); // should be P=20.00
        CPPUNIT_ASSERT( TW > 9.9 &&  TW < 10.1); // should be TW=10.000
        //10 deep
        YY = 10.000;
        CalculateSectionFlowParameters(channel,YY,&A,&P,&TW);
        CPPUNIT_ASSERT( A > 99.9 &&  A < 100.1); // should be A=100.00
        CPPUNIT_ASSERT( P > 29.9 &&  P < 30.1 ); // should be P=30.00
        CPPUNIT_ASSERT( TW > 9.9 &&  TW < 10.1); // should be TW=10.000
      }
      void testTrapSectionParameters()
      {
        ChannelShape channel;
        channel.CrossSectionShape = TRAP;
        // NOTE FOR TRAPEZIUM WIDTH IS BOTTOM WIDTH
        channel.Width = 10.000;
        double YY = 5.000;
        double A = 0;
        double P = 0;
        double TW = 0;
        //5m deep (almost rectangle)
        channel.SideWallAngle = 89.999;
        CalculateSectionFlowParameters(channel,YY,&A,&P,&TW);
        CPPUNIT_ASSERT( A > 49.9 &&  A < 50.1); // should be A=50.00
        CPPUNIT_ASSERT( P > 19.9 &&  P < 20.1); // should be P=20.00
        CPPUNIT_ASSERT( TW > 9.9 &&  TW < 10.1); // should be TW=10.000
       
       //5m deep (45 degrees)
        channel.SideWallAngle = 45.0;
        CalculateSectionFlowParameters(channel,YY,&A,&P,&TW);
        CPPUNIT_ASSERT( A > 74.9 &&  A < 75.1); // should be A=75.00
        CPPUNIT_ASSERT( P > 24.0 &&  P < 24.2); // should be P=24.14
        CPPUNIT_ASSERT( TW > 19.9 &&  TW < 20.1); // should be TW=20.000
      }
};