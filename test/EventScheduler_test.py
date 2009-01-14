#!/usr/bin/env python

import unittest

import weakref

import _gfrd as mod


class Delegate( object ):

    def __init__( self, obj, method ):
        self.obj = weakref.proxy( obj )
        self.method = method

    def __call__( self, arg ):
        return self.method( self.obj, arg )



class TestEvent:

    def __init__( self ):
        self.fired = 0
        self.dt = 1.0

    def fire( self, arg ):
        self.fired += 1
        return self.dt


class TestEvent2:

    def __init__( self, scheduler ):
        self.fired = 0
        self.dt = 1.0
        self.scheduler = scheduler

    def fire( self, arg ):
        self.fired += 1
        self.newevent = TestEvent()
        self.scheduler.addEvent( self.scheduler.getTime(), 
                                 Delegate( self.newevent, TestEvent.fire ), 
                                 () )
        return self.dt

class TestEvent3:

    def __init__( self, scheduler, peer ):
        self.fired = 0
        self.dt = 1.0
        self.scheduler = scheduler
        self.peer = peer

    def fire( self, arg ):
        self.fired += 1
        self.scheduler.updateEventTime( self.peer.id, self.scheduler.getTime() )
        return self.dt

class TestEvent4:

    def __init__( self, scheduler ):
        self.fired = 0
        self.dt = 1.0
        self.scheduler = scheduler

    def fire( self, arg ):
        self.fired += 1
        self.event1 = TestEvent()
        self.event2 = TestEvent()
        self.scheduler.addEvent( self.scheduler.getTime(), 
                                 Delegate( self.event1, TestEvent.fire ), 
                                 () )
        self.scheduler.addEvent( self.scheduler.getTime(), 
                                 Delegate( self.event2, TestEvent.fire ), 
                                 () )
        return self.dt




    

class EventSchedulerTestCase( unittest.TestCase ):

    def setUp(self):
        pass

    def tearDown(self):
        pass
    
    def testInstantiation(self):
        scheduler = mod.EventScheduler()
        self.failIf( scheduler == None )

    def testEmptyState(self):
        scheduler = mod.EventScheduler()
        self.failIf( scheduler.getSize() != 0 )
        self.failIf( scheduler.getTime() != 0.0 )
        # what if getTopEvent() are called here?

    def testOneEvent(self):
        scheduler = mod.EventScheduler()

        event = TestEvent()
        id = scheduler.addEvent( 0.0, Delegate( event, TestEvent.fire ), () )
        self.failIf( scheduler.getTime() != 0.0 )
        self.failIf( scheduler.getTopTime() != 0.0 )
        self.failIf( scheduler.getTopID() != id )
        
        scheduler.step()
        self.failIf( scheduler.getTime() != 0.0 )
        self.failIf( scheduler.getTopTime() != 1.0 )
        self.failIf( scheduler.getTopTime() != 
                     scheduler.getTopEvent().getTime() )
        self.failIf( scheduler.getTopID() != id )

    def testTwoEventsSameTime1(self):
        scheduler = mod.EventScheduler()

        event1 = TestEvent()
        event2 = TestEvent()
        event1.dt = 1.0
        event2.dt = 0.5
        scheduler.addEvent( 0.0, Delegate( event1, TestEvent.fire ), () )
        id2 = scheduler.addEvent( 0.0, Delegate( event2, TestEvent.fire ), () )
        self.assertEqual( 2, scheduler.getSize() )
        self.failIf( scheduler.getTime() != 0.0 )
        self.failIf( scheduler.getTopTime() != 0.0 )
        
        scheduler.step()
        self.failIf( scheduler.getTime() != 0.0 )
        self.failIf( scheduler.getTopTime() != 0.0 )

        scheduler.step()
        self.failIf( scheduler.getTopTime() != 0.5 )
        self.failIf( scheduler.getTopTime() != 
                     scheduler.getTopEvent().getTime() )
        self.failIf( scheduler.getTopID() != id2 )


    def testTwoEventsSameTime2(self):
        scheduler = mod.EventScheduler()

        event1 = TestEvent4( scheduler )
        event1.dt = 0.5
        id1 = scheduler.addEvent( 0.0, Delegate( event1, TestEvent4.fire ), () )
        self.assertEqual( 1, scheduler.getSize() )
        self.failIf( scheduler.getTime() != 0.0 )
        self.failIf( scheduler.getTopTime() != 0.0 )
        
        scheduler.step()
        self.assertEqual( 3, scheduler.getSize() )
        self.failIf( scheduler.getTime() != 0.0 )
        self.failIf( scheduler.getTopTime() != 0.0 )

        scheduler.step()
        scheduler.step()
        self.assertEqual( 0.5, scheduler.getTopTime() )
        self.assertEqual( scheduler.getTopID(), id1 )



    def testEventPop(self):
        scheduler = mod.EventScheduler()

        event1 = TestEvent()
        event1.dt = - 1.0

        scheduler.addEvent( 0.0, Delegate( event1, TestEvent.fire ), () )
        self.failIf( scheduler.getSize() != 1 )
        
        scheduler.step()

        self.failIf( scheduler.getSize() != 0 )
        self.failIf( scheduler.getTime() != 0.0 )

        
    def testEventPop2(self):

        scheduler = mod.EventScheduler()

        event1 = TestEvent()
        event2 = TestEvent()
        event1.dt = 1.0
        event2.dt = -1.0

        id1 = scheduler.addEvent( 0.0, Delegate( event1, TestEvent.fire ), () )
        scheduler.addEvent( 0.5, Delegate( event2, TestEvent.fire ), () )
        self.failIf( scheduler.getSize() != 2 )
        
        scheduler.step()
        self.failIf( scheduler.getTime() != 0.0 )
        scheduler.step()

        self.failIf( scheduler.getTime() != 0.5 )
        self.failIf( scheduler.getSize() != 1 )
        self.failIf( scheduler.getTopID() != id1 )
        self.failIf( scheduler.getTopTime() != 1.0 )


    def testEventCreationDuringStepping(self):

        scheduler = mod.EventScheduler()

        event1 = TestEvent()
        event2 = TestEvent2( scheduler )
        event1.dt = 1.0
        event2.dt = 1.0

        id1 = scheduler.addEvent( 0.0, Delegate( event1, TestEvent.fire ), () )
        id2 = scheduler.addEvent( 0.5, Delegate( event2, TestEvent2.fire ), () )
        self.failIf( scheduler.getSize() != 2 )
        
        scheduler.step()
        self.failIf( scheduler.getTime() != 0.0 )
        scheduler.step()

        self.failIf( scheduler.getTime() != 0.5 )
        self.failIf( scheduler.getSize() != 3 )
        self.failIf( scheduler.getTopID() == id1 or
                     scheduler.getTopID() == id2 )
        self.failIf( scheduler.getTopTime() != 0.5 )


    def testEventCreationDuringSteppingWithPop(self):

        scheduler = mod.EventScheduler()

        event1 = TestEvent()
        event2 = TestEvent2( scheduler )
        event1.dt = 1.0
        event2.dt = - 1.0

        id1 = scheduler.addEvent( 0.0, Delegate( event1, TestEvent.fire ), () )
        id2 = scheduler.addEvent( 0.5, Delegate( event2, TestEvent2.fire ), () )
        self.failIf( scheduler.getSize() != 2 )
        
        scheduler.step()
        self.failIf( scheduler.getTime() != 0.0 )
        scheduler.step()

        self.failIf( scheduler.getTime() != 0.5 )
        self.failIf( scheduler.getSize() != 2 )
        self.failIf( scheduler.getTopID() == id1 or
                     scheduler.getTopID() == id2 )
        self.failIf( scheduler.getTopTime() != 0.5 )


    def testEventUpdateDuringStepping(self):

        scheduler = mod.EventScheduler()

        event1 = TestEvent()
        event2 = TestEvent3( scheduler, event1 ) 
        event1.dt = 1.0
        event2.dt = 1.0

        event1.id = scheduler.addEvent( 0.0, Delegate( event1, TestEvent.fire ), () )

        id2 = scheduler.addEvent( 0.5, Delegate( event2, TestEvent3.fire ), () )
        self.assertEqual( 2, scheduler.getSize() )
        
        self.failIf( scheduler.getTopID() != event1.id )
        scheduler.step()
        self.assertEqual( 0.0, scheduler.getTime() )

        # event2 updates event1 to step immediately
        self.failIf( scheduler.getTopID() != id2 )
        scheduler.step()

        self.assertEqual( 0.5, scheduler.getTime() )
        self.assertEqual( 2, scheduler.getSize() )
        self.failIf( scheduler.getTopID() != event1.id )
        self.failIf( scheduler.getTopTime() != 0.5 )




if __name__ == "__main__":
    unittest.main()
