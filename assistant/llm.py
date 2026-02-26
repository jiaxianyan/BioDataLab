import os
from typing import Literal, Optional

import openai
from langchain_anthropic import ChatAnthropic

# from langchain_aws import ChatBedrock
from langchain_core.language_models.chat_models import BaseChatModel
from langchain_google_genai import (
    ChatGoogleGenerativeAI,
    HarmBlockThreshold,
    HarmCategory,
)
from langchain_ollama import ChatOllama
from langchain_openai import AzureChatOpenAI, ChatOpenAI
from langchain_xai import ChatXAI
from langchain_deepseek import ChatDeepSeek

SourceType = Literal["OpenAI", "AzureOpenAI", "Anthropic", "Ollama", "Gemini", "Bedrock", "Groq", "Custom"]
ALLOWED_SOURCES: set[str] = set(SourceType.__args__)

API_KEY = ""
BASE_URL = ""

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

    if model.startswith('gpt') or model.startswith('qwen'):
        llm = ChatOpenAI(model=model, 
                        base_url=BASE_URL,
                        openai_api_key=API_KEY,
                        temperature=temperature
                        )
        
        llm_with_stop = llm.bind(stop=stop_sequences)
        return llm_with_stop

    elif 'kimi' in model or 'deepseek' in model:
        llm = ChatOpenAI(model=model, 
                        base_url=BASE_URL,
                        openai_api_key=API_KEY,
                        temperature=1
                        )
        
        llm_with_stop = llm.bind(stop=stop_sequences)
        return llm_with_stop    
    
    elif model.startswith('claude'):
        llm = ChatAnthropic(
            model=model,  
            anthropic_api_url=BASE_URL, 
            anthropic_api_key=API_KEY,  
            temperature=temperature
        )
        llm_with_stop = llm.bind(stop=stop_sequences)
        return llm_with_stop
    
    elif model.startswith('gemini'):
        llm = ChatGoogleGenerativeAI(
                                      model=model,
                                      transport="rest",  # 非常关键
                                      client_options= {"api_endpoint": BASE_URL},
                                      google_api_key=API_KEY,
                                      thinking_level="low",
                                      safety_settings={
                                                        HarmCategory.HARM_CATEGORY_DANGEROUS_CONTENT: HarmBlockThreshold.BLOCK_NONE,
                                                    }
                                     )
        llm_with_stop = llm.bind(stop=stop_sequences)
        return llm_with_stop
    else:
        llm = ChatOpenAI(model=model, 
                        base_url=BASE_URL,
                        openai_api_key=API_KEY,
                        temperature=temperature, 
                        stop_sequences=stop_sequences)
                        
        llm_with_stop = llm.bind(stop=stop_sequences)
        return llm_with_stop