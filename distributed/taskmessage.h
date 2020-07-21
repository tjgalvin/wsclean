#ifndef TASK_MESSAGE_H
#define TASK_MESSAGE_H

struct TaskMessage {
  enum Type { Finish, GriddingRequest, GriddingResult } type;
  size_t bodySize = 0;
};

#endif
