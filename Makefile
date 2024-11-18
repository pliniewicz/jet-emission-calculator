CC := gcc

SRC_DIR   := src
BIN_DIR := bin
EXE 	  := $(BIN_DIR)/jets
CPPFLAGS := -Iinclude -MMD -MP
CFLAGS   := -Wall     
LDFLAGS  := -Llib     
LDLIBS   := -lm -lgsl -lcuba

SRCS := $(shell find $(SRC_DIR) -name '*.c')
OBJS := $(subst $(SRC_DIR), $(BIN_DIR), $(SRCS:.c=.o))

all : $(OBJS) $(EXE)

$(EXE) : $(OBJS) | $(BIN_DIR)
	@echo "------ Make $(EXE) ------"
	rm -f $(EXE)
	gcc $(CFLAGS) $(CPPFLAGS) $(LDFLAGS) $(LDLIBS) -o $(EXE) $(OBJS)

$(BIN_DIR)/%.o : $(SRC_DIR)/%.c | $(BIN_DIR)
	@echo "------ Make $(@) ------"
	rm -f $@
	gcc $(CFLAGS) $(CPPFLAGS) $(LDFLAGS) $(LDLIBS)  -c -o $@ $<

$(BIN_DIR):
	mkdir -p $(BIN_DIR)

-include $(BIN_DIR)/*.d

clean:
	rm -rf $(BIN_DIR)/*
