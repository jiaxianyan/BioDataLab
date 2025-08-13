import os
from typing import Literal, Optional

import openai
from langchain_anthropic import ChatAnthropic

# from langchain_aws import ChatBedrock
from langchain_core.language_models.chat_models import BaseChatModel
from langchain_google_genai import ChatGoogleGenerativeAI
from langchain_ollama import ChatOllama
from langchain_openai import AzureChatOpenAI, ChatOpenAI

SourceType = Literal["OpenAI", "AzureOpenAI", "Anthropic", "Ollama", "Gemini", "Bedrock", "Groq", "Custom"]
ALLOWED_SOURCES: set[str] = set(SourceType.__args__)


def get_llm(
    model: str = "claude-3-5-sonnet-20241022",
    temperature: float = 0.7,
    stop_sequences: list[str] | None = None,
    source: SourceType | None = None,
    base_url: str | None = None,
    api_key: str = "EMPTY",
) -> BaseChatModel:
    """
    Get a language model instance based on the specified model name and source.
    This function supports models from OpenAI, Azure OpenAI, Anthropic, Ollama, Gemini, Bedrock, and custom model serving.
    Args:
        model (str): The model name to use
        temperature (float): Temperature setting for generation
        stop_sequences (list): Sequences that will stop generation
        source (str): Source provider: "OpenAI", "AzureOpenAI", "Anthropic", "Ollama", "Gemini", "Bedrock", or "Custom"
                      If None, will attempt to auto-detect from model name
        base_url (str): The base URL for custom model serving (e.g., "http://localhost:8000/v1"), default is None
        api_key (str): The API key for the custom llm
    """

    # Create appropriate model based on source
    return ChatOpenAI(model=model, 
                        base_url="https://aihubmix.com/v1",
                        openai_api_key="sk-dlrre7E4eIX60dKU1eFdDb54F82b42AeB28c04716d211651",
                        temperature=temperature, 
                        max_tokens=8192,
                        stop_sequences=stop_sequences)