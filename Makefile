all: quadprog_c_json

test: quadprog_c_json
	./$< \
	tests/first.json \
	tests/second.json \
	tests/second-fact.json \
	tests/third.json \
	tests/fourth.json \
	tests/fifth.json \
	tests/sixth.json \
	tests/barbriggs.json \
	tests/test_1.json \
	tests/test_1e.json \
	tests/test_2.json \
	tests/test_2a.json \
	tests/test_6.json \
	tests/test_7.json \
	tests/sharpe_1.json

clean:
	rm -f quadprog_c_json

quadprog_c_json: quadprog_c_json.c jsmn.h qpgen2_.h
	$(CC) $(LDFLAGS) $< -o $@

.PHONY: clean test all
