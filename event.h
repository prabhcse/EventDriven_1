#ifndef EV_H
#define EV_H
// Prabhmanmeet Singh
// Like a Linked List
struct Event {
  double time;            // Time at which Event takes place Arrival or Departure
  int type;               // Type of Event (Arrival, Departure)
  int job;			  	  // Type of Job (Admin - 1, User - 2, Dep - 3)
  Event* next;            // Points to next event in list
  Event(double t, int i, int j) {
    time = t;
    type = i;
    job = j;
    next = 0;
  }
};

class EventList {
  public:
  Event* head;           // Points to first Event in EventList
  int event_count;       // Total number of Events in EventList
  ~EventList() { clear();}
  EventList() { event_count = 0; head = 0;}
  void insert(double time, int type, int job);  // Insert new event into EventList
  Event* get();                        				 // Returns first Event in EventList
  void clear();                        				 // Removes all Events from EventList
  Event* remove(int type);             // Returns first Event of given type and job type 
};

#endif
