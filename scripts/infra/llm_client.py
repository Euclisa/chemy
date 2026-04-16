import threading
import os

from contextlib import contextmanager

try:
    from openai import OpenAI
except ImportError:
    OpenAI = None


class LLMClient:

    class CompletionTokensLimitReached(Exception):
        def __init__(self, message="Reached limit of max completion tokens"):
            super().__init__(message)

    def __init__(self, api_key=None, logger=None):
        if api_key is None:
            api_key = os.getenv("OPENROUTER_API_KEY")

        self.client = None
        if OpenAI is not None:
            self.client = OpenAI(base_url="https://openrouter.ai/api/v1", api_key=api_key)

        self.completion_tokens_total = 0
        self.input_tokens_total = 0
        self.tokens_total_lock = threading.Lock()
        self._max_ctt = None
        self.logger = logger

    @contextmanager
    def restrict_ctt(self, max_ctt):
        self._max_ctt = max_ctt
        try:
            yield
        finally:
            self._max_ctt = None

    def _create_messages_block(self, message):
        messages = [
            {"role": "system", "content": ""},
            {"role": "user", "content": message}
        ]
        return messages

    def _fetch_answer(self, messages, model, reasoning_effort="medium", max_completion_tokens=10000):
        if self.client is None:
            raise ModuleNotFoundError("openai is required to use LLMClient")

        if self._max_ctt is not None and self.completion_tokens_total >= self._max_ctt:
            raise self.CompletionTokensLimitReached()

        completion = self.client.chat.completions.create(
            model=model,
            messages=messages,
            reasoning_effort=reasoning_effort,
            max_completion_tokens=max_completion_tokens
        )

        with self.tokens_total_lock:
            self.input_tokens_total += completion.usage.prompt_tokens
            self.completion_tokens_total += completion.usage.completion_tokens

        return completion.choices[0].message.content

    def fetch_answer_str(self, message, model, reasoning_effort="medium", max_completion_tokens=10000):
        messages = self._create_messages_block(message)
        return self._fetch_answer(messages, model, reasoning_effort=reasoning_effort, max_completion_tokens=max_completion_tokens)

    def get_future_result(self, future, executor):
        try:
            result = future.result()
        except self.CompletionTokensLimitReached:
            executor.shutdown(wait=False)
            raise
        except Exception as e:
            if self.logger:
                self.logger.log_err(f"Exception occured: {e}")
            return None

        return result

    @contextmanager
    def track_usage(self):
        start_ct = self.completion_tokens_total
        start_it = self.input_tokens_total
        usage = {}
        try:
            yield usage
        finally:
            usage['completion_tokens'] = self.completion_tokens_total - start_ct
            usage['input_tokens'] = self.input_tokens_total - start_it
